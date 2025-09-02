from __future__ import annotations

import os
import sys
import contextlib
import glob
import warnings
import io
import math
import re
import json
from pathlib import Path
import functools
from typing import (
    Any,
    Sequence,
    List,
    Tuple,
    Dict,
    Optional,
    Union,
    TextIO,
    Iterator,
    Iterable,
    ContextManager,
    ClassVar,
)
from functools import cached_property
from dataclasses import dataclass
from multiprocessing.pool import ThreadPool

import pyjess

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from enzymm.utils import chunks, ranked_argsort, DummyPool
from enzymm.mcsa_info import load_mcsa_catalytic_residue_homologs_info
from enzymm.mcsa_info import ReferenceCatalyticResidue, NonReferenceCatalyticResidue


@dataclass(frozen=True)
class Vec3:
    """
    Class for storing 3D vectors in XYZ.

    Attributes:
        x : `float` X-Coordinate value
        y : `float` Y-Coordinate value
        z : `float` Z-Coordinate value

    """

    x: float
    y: float
    z: float

    def __post_init__(self):
        if math.isnan(self.x) or math.isnan(self.y) or math.isnan(self.z):
            raise ValueError(
                "Cannot create a Vec3 with NaN values. Likely the Jess superposition failed."
            )

    @classmethod
    def from_xyz(cls, item: Any) -> Vec3:
        """
        Create a `Vec3` instance from an object with x,y,z attributes.

        Args:
            item `any` : Any object with x, y, z attributes

        Returns:
            `Vec3` instance
        """
        return Vec3(item.x, item.y, item.z)

    @property
    def norm(self) -> float:
        """
        `float`: The vector norm (root of sum of squares)
        """
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self) -> Vec3:
        """
        `Vec3`: The vector devided by its norm
        """
        norm = self.norm
        if norm == 0:  # zero vector
            return Vec3.from_xyz(self)
        return self / norm

    def __matmul__(self, other: Vec3) -> float:
        """
        Overloads the @ operator to perform dot product between two `Vec3` vectors.

        Args:
            other: `Vec3` instance

        Returns:
            `Vec3`
        """
        if isinstance(other, Vec3):
            return self.x * other.x + self.y * other.y + self.z * other.z
        raise TypeError(f"Expected Vec3, got {type(other).__name__}")

    def __add__(self, other: Union[int, float, Vec3]) -> Vec3:
        """
        Overloads the + operator to add either an `int`, `float` or other `Vec3` to a `Vec3`.

        Args:
            other: `int` | `float` | `Vec3` argument to add

        Returns:
            `Vec3`
        """
        if isinstance(other, Vec3):
            return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x + other, self.y + other, self.z + other)
        raise TypeError(f"Expected int, float or Vec3, got {type(other).__name__}")

    def __truediv__(self, other: Union[int, float, Vec3]) -> Vec3:
        """
        Overloads the / operator to divide a `Vec3` object by either an `int`, `float` or other `Vec3`.

        Args:
            other: `int` | `float` | `Vec3` argument to divide by

        Returns:
            `Vec3`
        """
        if isinstance(other, Vec3):
            return Vec3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x / other, self.y / other, self.z / other)
        raise TypeError(f"Expected int, float or Vec3, got {type(other).__name__}")

    def __sub__(self, other: Union[int, float, Vec3]) -> Vec3:
        """
        Overloads the - operator to subtract either an `int`, `float` or other `Vec3` from a `Vec3`.

        Args:
            other: `int` | `float` | `Vec3` argument to subtract

        Returns:
            `Vec3`
        """
        if isinstance(other, Vec3):
            return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x - other, self.y - other, self.z - other)
        raise TypeError(f"Expected int, float or Vec3, got {type(other).__name__}")

    # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    def angle_to(self, other: Vec3) -> float:
        """
        Returns the angle in radians between two `Vec3` vectors.

        Args:
            other: `Vec3` to which to calculate the angle

        Returns:
            `float`: angle in radians

        """
        dot_product = self.normalize() @ other.normalize()
        if -1 <= dot_product <= 1:
            return math.acos(dot_product)
        else:
            # due to numerical errors two identical vectors may have a dot_product not exactly 1
            if math.isclose(dot_product, 1, rel_tol=1e-5):
                return 0
            # same but with opposite vectors
            elif math.isclose(dot_product, -1, rel_tol=1e-5):
                return math.pi
            else:
                raise ValueError(
                    f"ArcCos is not defined outside [-1,1]. self.vec is {[self.x, self.y, self.z]}, other vec is {[other.x, other.y, other.z]}"
                )


@dataclass(frozen=True, init=False)
class Residue:
    """
    Class for storing template residues (defined as 3 atoms) with relevant information.
    """

    _atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]
    _vec: Vec3
    _indices: Tuple[int, int]

    def __init__(
        self,
        atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom],
    ):
        """
        Initilaize a `Residue` instance from a triplet of `~pyjess.TemplateAtom` objects.

        Args:
            atoms: `tuple` of 3 `~pyjess.TemplateAtom` instances

        Returns:
            `Residue`
        """

        vec, indices = self.calc_residue_orientation(atoms)
        object.__setattr__(self, "_atoms", atoms)
        object.__setattr__(self, "_vec", vec)
        object.__setattr__(self, "_indices", indices)

    @staticmethod
    def calc_residue_orientation(
        atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]
    ) -> Tuple[Vec3, Tuple[int, int]]:
        """
        Method to calculate the residue orientation depending on the residue type.

        Args:
            atoms: `tuple` of 3 `~pyjess.TemplateAtom`

        Note:
            For symmetric atom triplets, the angle is calculated from the central atom to the midpoint between the two identical atom types.
            For non-symetric atom triplets, the angle is calculated between two atoms following the axis of polarization.

        Returns:
            `tuple`:  of `Residue.orientation_vector` and `Residue.orientation_vector_indices`
        """

        # dictionary in which the vectors from start to finish are defined for each aminoacid type
        # the orientation vector is calculated differently for different aminoacid types
        vector_atom_type_dict = {
            "GLY": ("C", "O"),
            "ALA": ("CA", "CB"),
            "VAL": ("CA", "CB"),
            "LEU": ("CA", "CB"),
            "ILE": ("CA", "CB"),
            "MET": ("CG", "SD"),
            "PHE": ("CZ", "mid"),
            "TYR": ("CZ", "OH"),
            "TRP": ("CZ2", "NE1"),
            "CYS": ("CB", "SG"),
            "PRO": ("C", "O"),
            "SER": ("CB", "OG"),
            "THR": ("CB", "OG1"),
            "ASN": ("CG", "OD1"),
            "GLN": ("CD", "OE1"),
            "LYS": ("CE", "NZ"),
            "ARG": ("CZ", "mid"),
            "HIS": ("CG", "ND1"),
            "ASP": ("CG", "mid"),
            "GLU": ("CD", "mid"),
            "PTM": ("CA", "CB"),
            "ANY": ("C", "O"),
        }
        try:
            vectup = vector_atom_type_dict[atoms[0].residue_names[0]]
        except KeyError as exc:
            raise KeyError(
                f"Residue orientation is not defined for the residue type {atoms[0].residue_names[0]}"
            ) from exc

        # In residues with two identical atoms, the vector is calculated from the middle atom to the mid point between the identical pair
        if vectup[1] == "mid":
            try:
                middle_index, middle_atom = next(
                    (index, atom)
                    for index, atom in enumerate(atoms)
                    if atom.atom_names[0] == vectup[0]
                )
                side1, side2 = [atom for atom in atoms if atom != middle_atom]
                midpoint = (Vec3.from_xyz(side1) + Vec3.from_xyz(side2)) / 2
                return midpoint - Vec3.from_xyz(middle_atom), (middle_index, 9)
            except StopIteration:
                raise ValueError(
                    f"Failed to find middle atom for amino-acid {atoms[0].residue_names[0]!r}"
                ) from None

        else:  # from first atom to second atom
            try:
                first_atom_index, first_atom = next(
                    (index, atom)
                    for index, atom in enumerate(atoms)
                    if atom.atom_names[0] == vectup[0]
                )
            except StopIteration:
                raise ValueError(
                    f"Failed to find first atom for amino-acid {atoms[0].residue_names[0]!r}"
                ) from None
            try:
                second_atom_index, second_atom = next(
                    (index, atom)
                    for index, atom in enumerate(atoms)
                    if atom.atom_names[0] == vectup[1]
                )
            except StopIteration:
                raise ValueError(
                    f"Failed to find second atom for amino-acid {atoms[0].residue_names[0]!r}"
                ) from None
            return Vec3.from_xyz(second_atom) - Vec3.from_xyz(first_atom), (
                first_atom_index,
                second_atom_index,
            )

    @classmethod
    def construct_residues_from_atoms(
        cls, atoms: Iterable[pyjess.TemplateAtom]
    ) -> List[Residue]:

        residues = []
        for atom_triplet in chunks(atoms, 3):  # yield chunks of 3 atom lines each
            if len(atom_triplet) != 3:
                raise ValueError(
                    f"Failed to construct residues. Got only {len(atom_triplet)} ATOM lines"
                )
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = {
                (atom.residue_names[0], atom.chain_id, atom.residue_number)
                for atom in atom_triplet
            }
            if len(unique_residues) != 1:
                raise ValueError(
                    f"Failed to construct residues. Atoms of different match_mode, chains, residue types or residue numbers {unique_residues} found in Atom triplet"
                )
            residues.append(Residue(atom_triplet))

        return residues

    @property
    def atoms(
        self,
    ) -> Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]:
        """
        Get the `tuple` of three `~pyjess.TemplateAtom` describing the residue.
        """
        return self._atoms

    @property
    def residue_name(self) -> str:
        """
        `str`: Get the amino-acid type as three letter code from the first atom.
        """
        return self.atoms[0].residue_names[0]

    @property
    def allowed_residues(self) -> str:
        """
        `str`: Get the allowed residue types as string of single letter codes.
        """
        convert_to_single = {
            "ALA": "A",
            "CYS": "C",
            "ASP": "D",
            "GLU": "E",
            "PHE": "F",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LYS": "K",
            "LEU": "L",
            "MET": "M",
            "ASN": "N",
            "PRO": "P",
            "GLN": "Q",
            "ARG": "R",
            "SER": "S",
            "THR": "T",
            "VAL": "V",
            "TRP": "W",
            "TYR": "Y",
            "XXX": "X",
        }
        return "".join(set(convert_to_single[i] for i in self.atoms[0].residue_names))

    @property
    def specific(self) -> int:
        """
        `bool`: Atoms with a match_mode greater 100 are unspecific. <100 is specific.
        """
        return all([atm.match_mode < 100 for atm in self.atoms])

    @property
    def backbone(self) -> bool:
        """
        `bool`: True if an atom may match a backbone atom.

        True if the atom may match backbone atoms.
        Check if the atom has 'ANY' or 'XXX' in its residue_names attribute
        """

        # TODO check what X or XXX means in pyjess!
        return (
            "ANY" in self.atoms[0].residue_names or "XXX" in self.atoms[0].residue_names
        )

    @property
    def residue_number(self) -> int:
        """
        `int`: Get the pdb residue number from the first atom.
        """
        return self.atoms[0].residue_number

    @property
    def chain_id(self) -> str:
        """
        `str`: Get the pdb chain_id from the first atom.
        """
        return self.atoms[0].chain_id

    @property
    def orientation_vector(self) -> Vec3:
        """
        `Vec3`: Calculate the residue orientation vector according to the residue type.
        """
        return self._vec

    @property
    def orientation_vector_indices(self) -> Tuple[int, int]:
        """
        `tuple`: Return the indices of the atoms
        between which the orientation vector was calculated according to the residue type.
        """
        return self._indices


@dataclass(frozen=True)
class AnnotatedResidue(Residue):
    """
    Child class inheriting from `Residue` for M-CSA annotated template residues.

    Attributes:
        is_mutated : `bool` : Wether the residue in the M-CSA reference PDB structure was mutated.
        is_metal_ligand : `str` : Wether the residue in the M-CSA reference PDB structure was metal coordinating.
        roles : `Tuple[str]` : Tuple of EMO codes describing catalytic roles.
        reference_idx: `int` : Residue indentifier of the reference residue.
        has_ptm: `bool` : Wether the residue in the M-CSA reference PDB structure was post-translationally modified.

    """

    is_mutated: bool
    is_metal_ligand: bool
    roles: Tuple[str, ...]
    reference_idx: int
    has_ptm: bool


@dataclass(frozen=True)
class Cluster:
    """
    Class for storing template cluster information.

    Attributes:
        id: `int` Index of the template cluster
        member: `int` Member index within the template cluster
        size: `int` Total number of members in the template cluster
    """

    id: int
    member: int
    size: int

    def __post_init__(self):
        if self.member > self.size:
            raise ValueError("Cluster member cannot be greater than cluster size")


class Template(pyjess.Template):
    """
    Class for storing templates and associated information.

    Inherits and extends from `pyjess.Template`
    """

    _CATH_MAPPING: ClassVar[Dict[str, List[str]]]
    _EC_MAPPING: ClassVar[Dict[str, str]]
    _PDB_SIFTS: ClassVar[Dict[str, Dict[str, List[str]] | str]]

    def __init__(
        self,
        *,
        residues: Sequence[Residue],
        pdb_id: Optional[str] = None,
        mcsa_id: Optional[int] = None,
        id: Optional[str] = None,
        template_id_string: Optional[str] = None,
        cluster: Optional[Cluster] = None,
        uniprot_id: Optional[str] = None,
        organism: Optional[str] = None,
        organism_id: Optional[str] = None,
        resolution: Optional[float] = None,
        experimental_method: Optional[str] = None,
        enzyme_discription: Optional[str] = None,
        represented_sites: Optional[int] = None,
        ec: Iterable[str] = (),
        cath: Iterable[str] = (),
    ):
        """
        Initialize a template.

        Keyword Arguments:
            residues: `sequence` of `~Residue` instances
            id: `str` Internal Template ID string. Default `None`
            pdb_id: `str` The PDB ID of the template
            template_id_string: `str` String in the ID line of a template
            mcsa_id: `int` The M-CSA entry index form which the template was generated
            cluster: `Cluster` Instance of the template
            uniprot_id: `str` UniProt Identifier of the Protein from which the template was generated
            organism: `str` Organism name of the Protein from which the template was generated
            organism_id: `str` Taxonomic Identifier of the Organism of the Protein from which the template was generated
            resolution: `float` Resolution of the Protein Structure from which the template was generated
            experimental_method: `str` Experimental method by which the Protein Structure of the template was resolved
            enzyme_description: `str` Text Discription of the Protein from which the template was generated
            represented_sites: `int` The number of Enzymes which this template is representative for
            ec: `list` of EC numbers associated with Enzymes this template represents
            cath: `list` of CATH numbers associated with Enzymes this template represents

        NOTE:
            It is recommended not to pass an id string.
            If id is `None`, the id string will be set to:

            > {`~Template.effective_size`}-Residues_{~`Template.template_id_string`}_Cluster_{`Cluster.id`}-{`Cluster.member`}-{`Cluster.size`}

            This identifier string should be unique.

        NOTE:
            residues can be constructed from a list of ~pyjess.TemplateAtoms via the
            staticmethod ~Residue.construct_residues_from_atoms(atoms=atoms)

        Returns:
            `Template`
        """
        if len(residues) == 0:
            raise ValueError(
                "Tried creating an `Template` from an empty list of residues!"
            )

        for residue in residues:
            if not isinstance(residue, type(residues[0])):
                raise ValueError(
                    f"Tried creating a `Template` from different types of `Residue` objects. Got {type(residue) and type(residues[0])}"
                )

        self.residues = tuple(residues)
        super().__init__([atom for r in self.residues for atom in r._atoms], id=id)
        self.pdb_id = pdb_id
        self.mcsa_id = mcsa_id
        self.template_id_string = template_id_string
        self.cluster = cluster
        self.uniprot_id = uniprot_id
        self.organism = organism
        self.organism_id = organism_id
        self.resolution = resolution
        self.experimental_method = experimental_method
        self.enzyme_discription = enzyme_discription
        self.represented_sites = represented_sites
        self.enzyme_discription = enzyme_discription
        self.ec = tuple(sorted({*ec, *self._add_ec_annotations()}))
        self.cath = tuple(sorted({*cath, *self._add_cath_annotations()}))

    def _state(self) -> Tuple:
        """Used only for computing a hash and for equality comparisons"""
        return (
            tuple(self.residues),
            self.id,
            self.pdb_id,
            self.mcsa_id,
            self.template_id_string,
            self.cluster,
            self.effective_size,
            self.uniprot_id,
            self.organism,
            self.organism_id,
            self.resolution,
            self.experimental_method,
            self.enzyme_discription,
            self.represented_sites,
            self.ec,
            self.cath,
        )

    def __reduce_ex__(self, protocol):

        return (
            functools.partial(
                Template,
                residues=self.residues,
                pdb_id=self.pdb_id,
                mcsa_id=self.mcsa_id,
                id=self.id,
                template_id_string=self.template_id_string,
                cluster=self.cluster,
                uniprot_id=self.uniprot_id,
                organism=self.organism,
                organism_id=self.organism_id,
                resolution=self.resolution,
                experimental_method=self.experimental_method,
                enzyme_discription=self.enzyme_discription,
                represented_sites=self.represented_sites,
                ec=self.ec,
                cath=self.cath,
            ),
            (),
        )

    def __copy__(self) -> Template:
        return self.copy()

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Template):
            self_state = self._state()
            other_state = other._state()
            return self_state == other_state
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, Template):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self), self._state()))

    @classmethod  # reading from text
    def loads(
        cls,
        text: str,
        id: str | None = None,
        warn: bool = False,
    ) -> Template:
        """
        Load Template from `str`. Calls `Template.load()`

        Arguments:
            text: `str` of Template to load
            id: `str` or `None` Internal pyjess string which will superseed the ID string parsed from the template file. Default `None`
            warn: `bool` If warnings should be printed. Default `False`

        Returns:
            `Template`
        """
        return cls.load(io.StringIO(text), id=id, warn=warn)

    @classmethod  # reading from TextIO
    def load(
        cls,
        file: TextIO | Iterator[str] | str | os.PathLike[str],
        id: str | None = None,
        warn: bool = False,
    ) -> Template:
        """
        Overloaded load to parse a `pyjess.Template` and its associated info into an `Template` object

        Arguments:
            file: `file-like` object or `str` or `path-like` from which to load
            id: `str` or `None` Internal pyjess string which will superseed the ID string parsed from the template file. Default `None`
            warn: `bool` If warnings should be printed. Default `False`

        Returns:
            `Template`
        """
        atoms = []
        metadata: dict[str, object] = {"ec": list(), "cath": list()}

        _PARSERS = {
            "ID": cls._parse_template_id_string,
            "PDB_ID": cls._parse_pdb_id,
            "UNIPROT_ID": cls._parse_uniprot_id,
            "MCSA_ID": cls._parse_mcsa_id,
            "CLUSTER": cls._parse_cluster,
            "ORGANISM_NAME": cls._parse_organism_name,
            "ORGANISM_ID": cls._parse_organism_id,
            "RESOLUTION": cls._parse_resolution,
            "EXPERIMENTAL_METHOD": cls._parse_experimental_method,
            "EC": cls._parse_ec,
            "CATH": cls._parse_cath,
            "ENZYME": cls._parse_enzyme_discription,
            "REPRESENTING": cls._parse_represented_sites,
        }

        # NOTE
        # currently mutliple templates per file are not supported by this parser function.

        # NOTE
        # this parser does not check for the existance of a header like REMARK TEMPLATE

        try:
            context: ContextManager[TextIO] = open(os.fspath(file))  # type: ignore
        except TypeError:
            context = contextlib.nullcontext(file)  # type: ignore

        with context as f:
            seen_lines = set()
            for line in filter(str.strip, f):
                tokens = line.split()
                if tokens[0] == "REMARK":
                    if len(tokens) == 1:
                        continue  # skip lines which just REMARK and nothing else
                    parser = _PARSERS.get(tokens[1])
                    if parser is not None:
                        if len(tokens) < 3:
                            raise IndexError(
                                f"Expected some annotation after the REMARK {tokens[1]} flag"
                            )
                        elif tokens[2].upper() in ["NONE", "?", "NAN", "NA"]:
                            continue
                        parser(tokens, metadata, warn=warn)
                elif tokens[0] == "ATOM":
                    if line in seen_lines:
                        raise ValueError("Duplicate Atom lines passed!")
                    seen_lines.add(line)
                    atoms.append(pyjess.TemplateAtom.loads(line))
                elif tokens[0] == "HETATM":
                    raise ValueError(
                        "Supplied template with HETATM record. HETATMs cannot be searched by Jess"
                    )
                else:
                    continue

        residues = Residue.construct_residues_from_atoms(atoms=atoms)

        return Template(
            residues=residues,
            id=id,
            **metadata,  # type: ignore # unpack everything parsed into metadata
        )

    def copy(self) -> Template:
        return type(self)(
            residues=self.residues,
            id=self.id,
            pdb_id=self.pdb_id,
            template_id_string=self.template_id_string,
            mcsa_id=self.mcsa_id,
            cluster=self.cluster,
            uniprot_id=self.uniprot_id,
            organism=self.organism,
            organism_id=self.organism_id,
            resolution=self.resolution,
            experimental_method=self.experimental_method,
            enzyme_discription=self.enzyme_discription,
            represented_sites=self.represented_sites,
            ec=self.ec,
            cath=self.cath,
        )

    def dumps(self) -> str:
        # """
        # Dump `AnnotatedTemplate` to a `str`. Calls `AnnotatedTemplate.dump()`
        # """
        buffer = io.StringIO()
        self.dump(buffer)
        return (
            buffer.getvalue()
        )  # returns entire content temporary file object as a string

    def dump(self, file: TextIO):
        # """
        # Dump `AnnotatedTemplate` to `file-like` object.

        # Arguments:
        #     file: `file-like` object to write to
        # """

        raise NotImplementedError

        file.write("REMARK TEMPLATE\n")
        if self.template_id_string:
            file.write(f"REMARK ID {self.template_id_string}\n")
        if self.pdb_id:
            file.write(f"REMARK PDB_ID {self.pdb_id}\n")
        if self.mcsa_id:
            file.write(f"REMARK MCSA_ID {self.mcsa_id}\n")
        if self.uniprot_id:
            file.write(f"REMARK UNIPROT_ID {self.uniprot_id}\n")
        if self.cluster:
            file.write(
                f"REMARK CLUSTER {'_'.join([str(self.cluster.id), str(self.cluster.member), str(self.cluster.size)])}\n"
            )
        if self.organism:
            file.write(f"REMARK ORGANISM_NAME {self.organism}\n")
        if self.organism_id:
            file.write(f"REMARK ORGANISM_ID {self.organism_id}\n")
        if self.resolution:
            file.write(f"REMARK RESOLUTION {self.resolution}\n")
        if self.experimental_method:
            file.write(f"REMARK EXPERIMENTAL_METHOD {self.experimental_method}\n")
        if self.enzyme_discription:
            file.write(f"REMARK ENZYME {self.enzyme_discription}\n")
        if self.represented_sites:
            file.write(
                f"REMARK REPRESENTING {self.represented_sites} CATALYTIC SITES\n"
            )
        if self.experimental_method:
            file.write(f"REMARK EXPERIMENTAL_METHOD {self.experimental_method}\n")
        if self.ec:
            file.write(f"REMARK EC {','.join(self.ec)}\n")
        if self.cath:
            file.write(f"REMARK CATH {','.join(self.cath)}")
        # file.write(f"REMARK SIZE {self.effective_size}\n")
        # file.write(f"REMARK DIMENSION {self.dimension}\n")
        # file.write(f"REMARK MULTIMERIC {self.multimeric}\n")
        # file.write(f"REMARK RELATIVE_ORDER {self.relative_order}\n")

        # for residue in self.residues:
        #     file.write(
        #         f"REMARK ORIENTATION_VECTOR OF RESIDUE {residue.residue_number}: between atom {residue.orientation_vector_indices[0]} and {residue.orientation_vector_indices[1]} {residue.orientation_vector.x:.3f} {residue.orientation_vector.y:.3f} {residue.orientation_vector.z:.3f}\n"
        #     )

        for residue in self.residues:
            for atom in residue.atoms:
                raise NotImplementedError(
                    "Still need to write a dump method for templates in PyJess"
                )
                # TODO atom.dump(file)

        file.write("END\n")

    @cached_property
    def effective_size(self) -> int:
        """`int`: The number of unique residues in the template, excluding backbone residues and unspecific residues."""
        # Number of residues as evaluated, the effective size
        # Effective size of a template is not necessarily equal to the number of atom triplets in a template:
        # Residues matching ANY amino acid type are not counted as they are too general
        # These have a value of 100 or hgiher in the second column indicating unspecific residues
        # Not all BB residues are unspecific! Some are targeted towards only Gly for example and thus have values < 100

        # Even if a Residue must match a particular AA, 6 atoms from a given residue may be selected
        # once for main chain and once for side chain
        # Therefore we only count unique template residues!

        # It seems that some templates even have the same 3-atom triplets at the same location twice. This I assume must be an error
        # again a reason to only count unique residues

        effective_size = 0
        for residue in self.residues:
            if (
                residue.specific and not residue.backbone
            ):  # type specific and not Backbone
                effective_size += 1
        return effective_size

    @cached_property
    def multimeric(
        self,
    ) -> bool:  # if the template is split across multiple protein chains
        """`bool`: True if the template contains residues from multiple protein chains."""
        return not all(
            res.chain_id == list(self.residues)[0].chain_id for res in self.residues
        )

        # NOTE
        # Some catalytic sites are at the interface between chains.
        # Sometimes these interaces are between heteromers, sometimes between homomers
        # Sometimes an enzyme is also only active as a multimer even though each chain
        # has its own catalytic site. These are allosteric effects.

        # Currently templates only describe momomeric catalytic sites
        # and catalytic sites at the interface between homomers.
        # Each template has its own UniProt ID

        # For some M-CSA entries the reference has an interfacial site
        # but some of the PDB homologs / Template derived structures are not interfacial
        # and vice versa. This means that assembly of the template structure
        # and the reference are not necessarily the same. Example: M-CSA 9 reference vs
        # PDB: 1o93 with residues 58 and 259
        # For this reason, we chose to ignore chain assignments in template matching.

    @cached_property
    def relative_order(
        self,
    ) -> List[int]:  # list with length of deduplicated template dimension
        """`list` of `int`: Relative order of residues in the template sorted by the pdb residue number.

        Note:
            This only works for non-multimeric templates. In this case returns '[0]'.
        """
        if self.multimeric:
            return [0]
        else:
            # Now extract relative template order
            return ranked_argsort([res.residue_number for res in self.residues])

    def _add_cath_annotations(self) -> List[str]:
        """`list`: Pull CATH Ids associated with that template from SIFTS and from the M-CSA"""
        cath_list = []
        if self.mcsa_id:
            cath_list.extend(
                self._CATH_MAPPING[str(self.mcsa_id)]
            )  # this is a bit inaccurate possibly ... shouldnt we pull via pdb_id....
        if self.pdb_id:
            pdbchains = set()
            for res in self.residues:
                pdbchains.add("{}{}".format(self.pdb_id, res.chain_id))

            # Iterating of all pdbchains which were part of the AnnotatedTemplate
            for pdbchain in pdbchains:
                # also include CATH annotations from PDB-SIFTS
                subdict = self._PDB_SIFTS.get(pdbchain)
                if subdict is not None:
                    sifts_caths = subdict.get("cath").copy()  # type: ignore
                    for cath in sifts_caths:
                        if cath != "?":
                            cath_list.append(cath)

        return cath_list

    def _add_ec_annotations(self) -> List[str]:
        """`list`: Pull EC Annotations associated with that template from SIFTS and from the M-CSA"""
        ec_list = []
        if self.mcsa_id is not None:
            ec_list.append(
                self._EC_MAPPING[str(self.mcsa_id)]
            )  # this is a bit inaccurate possibly ... shouldnt we pull via pdb_id....
        if self.pdb_id is not None:
            pdbchains = set()
            for res in self.residues:
                pdbchains.add("{}{}".format(self.pdb_id, res.chain_id))

            # Iterating of all pdbchains which were part of the AnnotatedTemplate
            for pdbchain in pdbchains:
                # also include EC annotations from PDB-SIFTS
                subdict = self._PDB_SIFTS.get(pdbchain)
                if subdict is not None:
                    sifts_ecs = subdict.get("ec").copy()  # type: ignore
                    for ec in sifts_ecs:
                        if ec != "?":
                            ec_list.append(ec)

        return ec_list

    @classmethod
    def _parse_pdb_id(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        if len(tokens[2]) != 4:
            raise ValueError(
                f"Found {tokens[2]} which has more than the expected 4 characters of a PDB_ID."
            )
        metadata["pdb_id"] = tokens[2].lower()

    @classmethod
    def _parse_template_id_string(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        metadata["template_id_string"] = tokens[2]

    @classmethod
    def _parse_uniprot_id(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        match = re.search(
            r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
            tokens[2],
        )
        if match:
            metadata["uniprot_id"] = match.group()
        else:
            raise ValueError(f"Did not find a valid UniProt ID, found {tokens[2]}")

    @classmethod
    def _parse_mcsa_id(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        try:
            metadata["mcsa_id"] = int(tokens[2])
        except ValueError as exc:
            raise ValueError(f"Did not find a M-CSA ID, found {tokens[2]}") from exc

    @classmethod
    def _parse_cluster(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        try:
            cluster = Cluster(*list(map(int, tokens[2].split("_"))))
            metadata["cluster"] = cluster
        except ValueError as exc:
            raise ValueError(
                f"Did not find a Cluster specification in the form <id>_<member>_<size>, found {tokens[2]}"
            ) from exc

    @classmethod
    def _parse_organism_name(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        metadata["organism"] = " ".join(tokens[2:])

    @classmethod
    def _parse_organism_id(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        metadata["organism_id"] = tokens[2]

    @classmethod
    def _parse_resolution(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        try:
            metadata["resolution"] = float(tokens[2])
        except ValueError as exc:
            raise ValueError(f"Ill-formatted pdb resolution: {tokens[2]}") from exc

    @classmethod
    def _parse_experimental_method(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        metadata["experimental_method"] = " ".join(tokens[2:])

    @classmethod
    def _parse_ec(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        matches = [
            match.group()
            for match in re.finditer(r"[1-7](\.(\-|\d{1,})){3}", tokens[2])
        ]
        non_cat_matches = [
            match.group()
            for match in re.finditer(r"[1-7](\.(\-|\d{1,}|n\d{1,})){3}", tokens[2])
        ]
        if matches:
            for ec in matches:
                if ec not in metadata["ec"]:  # type: ignore
                    metadata["ec"].append(ec)  # type: ignore
        elif non_cat_matches:
            for ec in non_cat_matches:
                if ec not in metadata["ec"]:  # type: ignore
                    metadata["ec"].append(ec)  # type: ignore
            if warn:
                warnings.warn(
                    f"Rare EC number(s) {[ec for ec in non_cat_matches]} presumed to be noncatalytic detected!"
                )
        else:
            raise ValueError(f"Did not find a valid EC number, found {tokens[2]}")

    @classmethod
    def _parse_cath(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        matches = [
            match.group()
            for match in re.finditer(r"[1-46](\.(\-|\d{1,})){3}", tokens[2])
        ]
        if matches:
            for cath in matches:
                if cath not in metadata["cath"]:  # type: ignore
                    metadata["cath"].append(cath)  # type: ignore
        else:
            raise ValueError(f"Did not find a valid CATH number, found {tokens[2]}")

    @classmethod
    def _parse_enzyme_discription(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        metadata["enzyme_discription"] = " ".join(tokens[2:])

    @classmethod
    def _parse_represented_sites(
        cls, tokens: List[str], metadata: dict[str, object], warn: bool = True
    ):
        try:
            metadata["represented_sites"] = int(tokens[2])
        except ValueError as exc:
            raise ValueError(
                f"Ill-formatted number of represented sites: {tokens[2]}"
            ) from exc


class AnnotatedTemplate(Template):

    residues: Tuple[AnnotatedResidue, ...]
    pdb_id: str
    mcsa_id: int

    def __init__(
        self,
        *,
        residues: Sequence[AnnotatedResidue],
        pdb_id: str,
        mcsa_id: int,
        number_of_mutated_residues: int,
        number_of_metal_ligands: Tuple[int, int],
        number_of_ptm_residues: Tuple[int, int],
        number_of_side_chain_residues: Tuple[int, int],
        total_reference_residues: int,
        assembly: int,
        id: Optional[str] = None,
        template_id_string: Optional[str] = None,
        cluster: Optional[Cluster] = None,
        uniprot_id: Optional[str] = None,
        organism: Optional[str] = None,
        organism_id: Optional[str] = None,
        resolution: Optional[float] = None,
        experimental_method: Optional[str] = None,
        enzyme_discription: Optional[str] = None,
        represented_sites: Optional[int] = None,
        ec: Iterable[str] = (),
        cath: Iterable[str] = (),
    ):
        """
        Initialize an annotated template with descriptions of catalytic activity from the M-CSA.

        Keyword Arguments:
            residues: `sequence` of `~Residue` instances
            id: `str` Internal Template ID string. Default `None`
            pdb_id: `str` The PDB ID of the template
            template_id_string: `str` String in the ID line of a template
            mcsa_id: `int` The M-CSA entry index form which the template was generated
            cluster: `Cluster` Instance of the template
            uniprot_id: `str` UniProt Identifier of the Protein from which the template was generated
            organism: `str` Organism name of the Protein from which the template was generated
            organism_id: `str` Taxonomic Identifier of the Organism of the Protein from which the template was generated
            resolution: `float` Resolution of the Protein Structure from which the template was generated
            experimental_method: `str` Experimental method by which the Protein Structure of the template was resolved
            enzyme_description: `str` Text Discription of the Protein from which the template was generated
            represented_sites: `int` The number of Enzymes which this template is representative for
            ec: `list` of EC numbers associated with Enzymes this template represents
            cath: `list` of CATH numbers associated with Enzymes this template represents
            number_of_mutated_residues: `int` The number of side chain specific residues which have been mutated relative to the reference
            number_of_metal_ligands: `tuple(`int`, `int`)` Number of metal chelating residues in the template and the reference
            number_of_ptm_residues: `tuple(`int`, `int`)` Number of post translationally modified residues in the template and the reference
            number_of_side_chain_residues: `tuple(`int`, `int`)` Number of side chain interacting residues in the template and the reference
            total_reference_residues: `int` Total number of residues (main and side chain) in the reference structure

        NOTE:
            In order for a template file to be loaded as an AnnotatedTemplate,
            it must need both an M-CSA id and a pdb_id.
            This pdb_id must be found in the PDB-homologs of the M-CSA!

        NOTE:
            It is recommended not to pass an id string.
            If id is `None`, the id string will be set to

            > {`~Template.effective_size`}-Residues_{~`Template.template_id_string`}_Cluster_{`Cluster.id`}-{`Cluster.member`}-{`Cluster.size`}

            This identifier string should be unique.

        Returns:
            `AnnotatedTemplate`
        """

        # if mcsa_id is None or pdb_id is None:
        #     raise ValueError("Missing mcsa_id or pdb_id for AnnotatedTemplate object.")

        super().__init__(
            residues=residues,
            pdb_id=pdb_id,
            mcsa_id=mcsa_id,
            id=id,
            template_id_string=template_id_string,
            cluster=cluster,
            uniprot_id=uniprot_id,
            organism=organism,
            organism_id=organism_id,
            resolution=resolution,
            experimental_method=experimental_method,
            enzyme_discription=enzyme_discription,
            represented_sites=represented_sites,
            ec=ec,
            cath=cath,
        )

        # Now we add in all the special information
        self.residues = tuple(residues)
        self.number_of_mutated_residues = number_of_mutated_residues
        self.number_of_metal_ligands = number_of_metal_ligands
        self.number_of_ptm_residues = number_of_ptm_residues
        self.number_of_side_chain_residues = number_of_side_chain_residues
        self.total_reference_residues = total_reference_residues
        self.assembly = assembly

    def __reduce_ex__(self, protocol):
        return (
            functools.partial(
                AnnotatedTemplate,
                residues=self.residues,
                pdb_id=self.pdb_id,
                mcsa_id=self.mcsa_id,
                number_of_mutated_residues=self.number_of_mutated_residues,
                number_of_metal_ligands=self.number_of_metal_ligands,
                number_of_ptm_residues=self.number_of_ptm_residues,
                number_of_side_chain_residues=self.number_of_side_chain_residues,
                total_reference_residues=self.total_reference_residues,
                assembly=self.assembly,
                id=self.id,
                template_id_string=self.template_id_string,
                cluster=self.cluster,
                uniprot_id=self.uniprot_id,
                organism=self.organism,
                organism_id=self.organism_id,
                resolution=self.resolution,
                experimental_method=self.experimental_method,
                enzyme_discription=self.enzyme_discription,
                represented_sites=self.represented_sites,
                ec=self.ec,
                cath=self.cath,
            ),
            (),
        )

    def _state(self) -> Tuple:

        residues = Residue.construct_residues_from_atoms(
            atoms=[atom for r in self.residues for atom in r._atoms]
        )
        # We use reuglar unanntoated residues to get the state!
        return (
            tuple(residues),
            self.id,
            self.pdb_id,
            self.mcsa_id,
            self.template_id_string,
            self.cluster,
            self.effective_size,
            self.uniprot_id,
            self.organism,
            self.organism_id,
            self.resolution,
            self.experimental_method,
            self.enzyme_discription,
            self.represented_sites,
            self.ec,
            self.cath,
        )

    def __copy__(self) -> AnnotatedTemplate:
        return self.copy()

    def __eq__(self, other: object) -> bool:
        if isinstance(other, AnnotatedTemplate):
            self_state = self._state()
            other_state = other._state()
            return self_state == other_state
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, AnnotatedTemplate):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self), self._state()))

    def copy(self) -> AnnotatedTemplate:
        return type(self)(
            residues=self.residues,
            id=self.id,
            pdb_id=self.pdb_id,
            template_id_string=self.template_id_string,
            mcsa_id=self.mcsa_id,
            cluster=self.cluster,
            uniprot_id=self.uniprot_id,
            organism=self.organism,
            organism_id=self.organism_id,
            resolution=self.resolution,
            experimental_method=self.experimental_method,
            enzyme_discription=self.enzyme_discription,
            represented_sites=self.represented_sites,
            ec=self.ec,
            cath=self.cath,
            number_of_mutated_residues=self.number_of_mutated_residues,
            number_of_metal_ligands=self.number_of_metal_ligands,
            number_of_ptm_residues=self.number_of_ptm_residues,
            number_of_side_chain_residues=self.number_of_side_chain_residues,
            total_reference_residues=self.total_reference_residues,
            assembly=self.assembly,
        )

    @classmethod  # reading from text
    def loads(
        cls,
        text: str,
        id: str | None = None,
        warn: bool = False,
        with_annotations: bool = True,
    ) -> Template | AnnotatedTemplate:
        """
        Load Template from `str`. Calls `Template.load()`

        Arguments:
            text: `str` of Template to load
            id: `str` or `None` Internal pyjess string which will superseed the ID string parsed from the template file. Default `None`
            warn: `bool` If warnings should be printed. Default `False`
            with_annotations: `bool` If True (default) M-CSA derived templates with a PDB-id and M-CSA id will be annotated with extra information.

        Returns:
            `Template|AnnotatedTemplate`
        """
        return cls.load(
            io.StringIO(text), id=id, warn=warn, with_annotations=with_annotations
        )

    @classmethod  # reading from TextIO
    def load(
        cls,
        file: TextIO | Iterator[str] | str | os.PathLike[str],
        id: str | None = None,
        warn: bool = False,
        with_annotations: bool = True,
    ) -> Template | AnnotatedTemplate:
        """
        Overloaded load to parse a `pyjess.Template` and its associated info into an `Template` object

        Arguments:
            file: `file-like` object or `str` or `path-like` from which to load
            id: `str` or `None` Internal pyjess string which will superseed the ID string parsed from the template file. Default `None`
            warn: `bool` If warnings should be printed. Default `False`
            with_annotations: `bool` If True (default) M-CSA derived templates with a PDB-id and M-CSA id will be annotated with extra information.

        Returns:
            `Template`|`AnnotatedTemplate`
        """

        template = Template.load(
            file=file,
            id=id,
            warn=warn,
        )

        if (
            with_annotations
            and template.pdb_id is not None
            and template.mcsa_id is not None
        ):
            annotated_residues, ann_dict = cls.derive_mcsa_annotations(template)

            return AnnotatedTemplate(
                residues=annotated_residues,
                id=template.id,
                pdb_id=template.pdb_id,
                mcsa_id=template.mcsa_id,
                template_id_string=template.template_id_string,
                cluster=template.cluster,
                uniprot_id=template.uniprot_id,
                organism=template.organism,
                organism_id=template.organism_id,
                resolution=template.resolution,
                experimental_method=template.experimental_method,
                enzyme_discription=template.enzyme_discription,
                represented_sites=template.represented_sites,
                ec=template.ec,
                cath=template.cath,
                number_of_mutated_residues=ann_dict["number_of_mutated_residues"],
                number_of_metal_ligands=ann_dict["number_of_metal_ligands"],
                number_of_ptm_residues=ann_dict["number_of_ptm_residues"],
                number_of_side_chain_residues=ann_dict["number_of_side_chain_residues"],
                total_reference_residues=ann_dict["total_reference_residues"],
                assembly=ann_dict["assembly_id"],
            )

        else:
            return template

    @staticmethod
    def derive_mcsa_annotations(
        template: Template,
    ) -> Tuple[List[AnnotatedResidue], Dict[str, Any]]:

        # NOTE
        # be careful how to interpret and handle this data.
        # Some templates are build from the reference pdb structure in the entry
        # Others, are built from a non-reference structure
        # Further, there are multimeric mcsa catalytic sites composed of:
        #   homo-mers (Example M-CSA 1, 4, ...)
        #   hetero-mers (Example M-CSA 5, 10, 11, ...)

        # NOTE
        # This only works if the template comes from the M-CSA (with id) and PDB id

        if template.mcsa_id is None or template.pdb_id is None:
            raise ValueError(
                "Tried annotating a template lacking either pdb_id or mcsa_id"
            )

        reference_homologs = set()
        annotated_residues = []
        assembly_ids = set()
        for residue in template.residues:
            ################### get reference residue ##################################
            try:
                template_pdbchain = CATALYTIC_RESIDUE_HOMOLOGS[template.mcsa_id][  # type: ignore
                    template.pdb_id + residue.chain_id  # type: ignore
                ]
            except KeyError:
                try:
                    template_pdbchain = CATALYTIC_RESIDUE_HOMOLOGS[template.mcsa_id][  # type: ignore
                        template.pdb_id + residue.chain_id[0]  # type: ignore
                    ]
                    # print(f"Only found a template pdb for {template.pdb_id} with pdbchain {residue.chain_id[0]} instead of {residue.chain_id}")
                    # these are all homo-mers except for 1olx
                    # but which is correctly assigned here too
                except KeyError:
                    raise KeyError(
                        f"Failed to find template pdb in catalytic residue homologs for M-CSA id {template.mcsa_id} and pdbchain {template.pdb_id+residue.chain_id}"  # type: ignore
                    ) from None

            assembly_ids.add(template_pdbchain.assembly)

            match_found = False
            for index, hom_residue in template_pdbchain.residues.items():
                # TODO check if it would even make a difference
                # TODO Ask ioannis about this: Should I check auth_resid or resid first?
                if residue.residue_number == hom_residue.auth_resid:
                    match_found = True
                    break
                elif residue.residue_number == hom_residue.resid:
                    match_found = True
                    break

            if not match_found:
                raise ValueError(
                    f"Missing a comparison residue for M-CSA id {template.mcsa_id} and pdbchain {template.pdb_id+residue.chain_id} for residue {residue.residue_name, residue.residue_number}"  # type: ignore
                ) from None

            # after the for loop breaks, index contains the residue index
            # and hom_residue contains the corresponding residue annotations
            if isinstance(hom_residue, ReferenceCatalyticResidue):
                reference_pdb = template_pdbchain
                ref_residue = hom_residue
            elif isinstance(hom_residue, NonReferenceCatalyticResidue):
                _, ref_pdbchain, _ = hom_residue.reference
                reference_pdb = CATALYTIC_RESIDUE_HOMOLOGS[template.mcsa_id][  # type: ignore
                    ref_pdbchain
                ]
                ref_residue = reference_pdb.residues[index]  # type: ignore

            reference_homologs.add(reference_pdb)

            ############################################################################

            # We check mutation only if
            # the type of the template residue matters
            # and it is interacting via its sidechain
            is_mutated = False
            if (
                residue.residue_name not in ["ANY", "PTM"]
                and residue.specific
                and ref_residue.function_location_abv is not None
            ):
                is_mutated = residue.residue_name != ref_residue.code

            # if template.mcsa_id == 661 and template.pdb_id == "4cyr" and ref_residue.resid == 51:
            #     print(ref_pdbchain, ref_residue.resid, ref_residue.code, ref_residue.ptm, bool(ref_residue.ptm))

            #     if isinstance(hom_residue, NonReferenceCatalyticResidue):
            #         print("is non_reference")

            annotated_residues.append(
                AnnotatedResidue(
                    _atoms=residue._atoms,
                    _vec=residue._vec,
                    _indices=residue._indices,
                    reference_idx=index,
                    is_metal_ligand="metal ligand" in ref_residue.roles_summary,
                    is_mutated=is_mutated,
                    has_ptm=bool(ref_residue.ptm),
                    roles=tuple(ref_residue.roles),
                )
            )

        if len(assembly_ids) == 1:
            assembly_id = list(assembly_ids)[0]
        else:
            raise ValueError(
                f"Got multiple assemblies for template from {template.pdb_id}"
            )

        ######################### Template level counts ################################
        # for templates where all the residues have annotations
        # we loop again over these residues to get template level counts
        number_side_chain_residues = 0
        number_metal_ligands = 0
        number_ptm = 0
        number_mutated = 0

        for residue in annotated_residues:
            if residue.specific:
                number_side_chain_residues += 1
            if residue.is_metal_ligand:
                number_metal_ligands += 1
            if residue.has_ptm:
                number_ptm += 1
            if residue.is_mutated:
                number_mutated += 1

        ################################################################################

        ############### Reference Homolog counts #######################################
        number_reference_residues = 0
        number_side_chain_reference_residues = 0
        number_metal_ligands_reference = 0
        number_ptm_reference = 0
        for reference_pdb in reference_homologs:
            for residue_id, residue in reference_pdb.residues.items():  # type: ignore
                if isinstance(residue, ReferenceCatalyticResidue):
                    # skip if the residue doesnt exist in the reference
                    if residue.code is None:
                        continue
                    number_reference_residues += 1
                    # function_location_abv is only set if it is NOT a side chain interaction
                    if not residue.function_location_abv:
                        number_side_chain_reference_residues += 1
                    if "metal ligand" in residue.roles_summary:
                        number_metal_ligands_reference += 1
                    if residue.ptm:
                        number_ptm_reference += 1

        ######################## assing Template properties/attributes #################

        return annotated_residues, {
            # number of mutated residues in the template versus the reference.
            "number_of_mutated_residues": number_mutated,
            # tuple(template, reference)
            "number_of_metal_ligands": (
                number_metal_ligands,
                number_metal_ligands_reference,
            ),
            # tuple(template, reference)
            "number_of_ptm_residues": (number_ptm, number_ptm_reference),
            # Side chain residues and ANY and PTM are discarded
            # tuple(template, reference)
            "number_of_side_chain_residues": (
                number_side_chain_residues,
                number_side_chain_reference_residues,
            ),
            # total number of residues in the reference pdb structure
            "total_reference_residues": number_reference_residues,
            "assembly_id": assembly_id,
        }


# Populate the mapping of MCSA IDs to CATH numbers so that it can be accessed
# by individual templates in the `Template.cath` property.
# Source: M-CSA which provides cath annotations for either residue homologs or for m-csa entries
with resource_files(__package__).joinpath("data/MCSA_CATH_mapping.json").open() as f:
    Template._CATH_MAPPING = json.load(f)

with resource_files(__package__).joinpath("data/MCSA_EC_mapping.json").open() as f:
    Template._EC_MAPPING = json.load(f)

# Source: CATH, EC and InterPro from PDB-SIFTS through mapping to the pdbchain
with resource_files(__package__).joinpath("data/pdb_sifts.json").open() as f:
    Template._PDB_SIFTS = json.load(f)

CATALYTIC_RESIDUE_HOMOLOGS = load_mcsa_catalytic_residue_homologs_info(
    resource_files(__package__).joinpath("data")  # type: ignore
)

# global MCSA_interpro_dict
# # dictonariy mapping M-CSA entries to Interpro Identifiers
# # Interpro Acceccesions at the Domain, Family and Superfamily level
# # are searched for the reference sequences of each M-CSA entry.
# # Note that an M-CSA entry may have multiple reference sequences

# TODO
# # add a list of cofactors associated with each EC number from Neeras List
# cofactors = set()
# if Template_EC in cofactor_dict:
#     cofactors.update(cofactor_dict[Template_EC])


# Find all template files which end in .pdb
def _get_paths_by_extension(directory_path: Path, extension: str) -> List[Path]:
    pattern = f"{directory_path}/**/*{extension}"
    file_paths = glob.glob(pattern, recursive=True)
    files = []
    for file in file_paths:
        files.append(Path(file))

    if files:
        return files
    else:
        raise FileNotFoundError(
            f"No template files with the {extension} extension found in the {directory_path.resolve()} directory"
        )


def load_templates(
    template_dir: Optional[Path] = None,
    warn: bool = False,
    verbose: bool = False,
    cpus: int = 1,  # faster with one thread
    # cpus: int = (
    #     len(os.sched_getaffinity(0)) if sys.platform == "linux" else os.cpu_count() or 1
    # ),
    with_annotations: bool = True,
) -> Iterator[Template | AnnotatedTemplate]:
    """
    Load templates from a given directory, recursively.

    Arguments:
        template_dir: `~pathlib.Path` | `None` Directory which to search recursively for files with the '.pdb' extension. By default, set to `None`, it will load templates included in this library.
        warn: `bool` If warnings about annoation issues in templates should be printed. Default `False`
        verbose: `bool` If loading should be verbose. Default `False`
        cpus: `int` The number of CPU threads to use. By default, use 1 thread. If <0, leave this number of threads free.
        with_annotations: `bool` If True (default) M-CSA derived templates with a PDB-id and M-CSA id will be annotated with extra information.

    Yields:
        `Template`|`AnnotatedTemplate`
    """
    if template_dir is None:
        template_dir = Path(
            str(resource_files(__package__).joinpath("jess_templates_20230210"))
        )

    elif isinstance(template_dir, Path):
        if not template_dir.is_dir():
            raise NotADirectoryError(
                f"The path {template_dir} doesnt exist or is not a directory!"
            )

    if verbose:
        print(f"Loading Template files from {str(template_dir.resolve())}")

    template_paths = _get_paths_by_extension(template_dir, ".pdb")

    def _load_and_annotate(template_path: Path, warn: bool, with_annotations: bool):
        try:
            with template_path.open() as f:
                return AnnotatedTemplate.load(
                    file=f, warn=warn, with_annotations=with_annotations
                )
        except ValueError as exc:
            raise ValueError(
                f"Passed Template file {template_path.resolve()} contained ATOM lines which are not in Jess Template format."
            ) from exc
        except KeyError as exc:
            raise ValueError(
                f"Passed Template file {template_path.resolve()} contained issues with some residues."
            ) from exc

    if cpus <= 0:
        os_cpu_count = (
            len(os.sched_getaffinity(0)) if sys.platform == "linux" else os.cpu_count()
        )
        if os_cpu_count is not None:
            cpus = max(1, os_cpu_count + cpus)
        else:
            cpus = 1

    pool: DummyPool | ThreadPool = DummyPool() if cpus == 1 else ThreadPool(cpus)
    with pool:
        args = [(path, warn, with_annotations) for path in template_paths]
        try:
            template_iterator = pool.starmap(_load_and_annotate, args)
            for loaded_template in template_iterator:
                yield loaded_template
        except ValueError as exc:
            raise exc


# TODO adapt this too for AnnotatedTemplate objects
def check_template(template: Template, warn: bool = True) -> bool:
    # TODO improve this and write tests
    if warn:
        checker = True

        # Raise warnings if some properties could not be annotated!
        if not template.ec:
            checker = False
            warnings.warn("Could not find EC number annotations")

        if not template.cath:
            checker = False
            warnings.warn("Could not find CATH annotations")

        if template.pdb_id:
            # check overlap between sifts mapping and CATH, EC annotations
            # Source: pdb to sifts mapping which maps CATH to pdb chain IDs and UniProt IDs, sifts also provides UniProt to EC mapping
            # Since Template may contain residues from multiple chains
            pdbchains = set()
            for res in template.residues:
                pdbchains.add(template.pdb_id + res.chain_id)

            # Iterating of all pdbchains which were part of the Template
            for pdbchain in pdbchains:
                # also include EC annotations from PDB-SIFTS
                subdict = template._PDB_SIFTS.get(pdbchain, None)
                if subdict:
                    # if template.cath and subdict['cath']:
                    #     if not set(template.cath) & set(subdict['cath']):
                    #         warnings.warn(f'Did not find an intersection of CATH domains as annotated by the M-CSA ID {template.mcsa_id} with {template.cath} versus PDF-SIFTS {template._PDB_SIFTS[pdbchain]['cath']} for template {filepath} from PDB ID {template.pdb_id}')
                    sifts_uniprot = subdict["uniprot_id"]  # type: ignore
                    if template.uniprot_id != sifts_uniprot:
                        checker = False
                        warnings.warn(
                            f"Different UniProt Accessions {template.uniprot_id} and {sifts_uniprot} found"
                        )

        return checker
    else:
        return True
