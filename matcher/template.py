from __future__ import annotations

import os
import glob
import warnings
import re
import json
from importlib.resources import files
from pathlib import Path
from typing import Any, Sequence, List, Tuple, Dict, Optional, Union, TextIO, Iterator, Iterable, ClassVar, IO
from functools import cached_property
from dataclasses import dataclass, field
import io
import math
from multiprocessing.pool import ThreadPool

import pyjess

from .utils import chunks, ranked_argsort, DummyPool

__all__ = [
    "Vec3",
    "Residue",
    "AnnotatedTemplate",
    "load_templates",
]

@dataclass(frozen=True)
class Vec3:
    '''Class for storing 3D vectors in XYZ'''
    x: float
    y: float
    z: float

    def __post_init__(self):
        if math.isnan(self.x) or math.isnan(self.y) or math.isnan(self.z):
            raise ValueError("Cannot create a Vec3 with NaN values. Likely the Jess superposition failed.")

    @classmethod
    def from_xyz(cls, item: Any) -> Vec3:
        return Vec3(item.x, item.y, item.z)

    @property
    def norm(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self) -> Vec3:
        norm = self.norm
        if norm == 0: # zero vector
            return Vec3.from_xyz(self)
        return self/norm

    def __matmul__(self, other: Vec3) -> float:
        """
        Overloads the @ operator to perform dot product between two vectors
        """
        if isinstance(other, Vec3):
            return self.x * other.x + self.y * other.y + self.z * other.z
        raise TypeError(f'Expected Vec3, got {type(other).__name__}')

    def __add__(self, other: Union[int, float, Vec3]) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x + other, self.y + other, self.z + other)
        raise TypeError(f'Expected int, float or Vec3, got {type(other).__name__}')

    def __truediv__(self, other: Union[int, float, Vec3]) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x / other, self.y / other, self.z / other)
        raise TypeError(f'Expected int, float or Vec3, got {type(other).__name__}')

    def __sub__(self, other: Union[int, float, Vec3]) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
        elif isinstance(other, (int, float)):
            return Vec3(self.x - other, self.y - other, self.z - other)
        raise TypeError(f'Expected int, float or Vec3, got {type(other).__name__}')

    # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    def angle_to(self, other: Vec3) -> float: 
        """ Returns the angle in radians between vectors 'v1' and 'v2':
        """
        dot_product = self.normalize() @ other.normalize()
        if -1 <= dot_product <= 1:
            return math.acos(dot_product)
        else:
            # due to numerical errors two identical vectors may have a dot_product not exactly 1
            if math.isclose(dot_product,1, rel_tol=1e-5): 
                return 0
            # same but with opposite vectors
            elif math.isclose(dot_product,-1, rel_tol=1e-5):
                return math.pi
            else:
                raise ValueError(f'ArcCos is not defined outside [-1,1]. self.vec is {[self.x, self.y, self.z]}, other vec is {[other.x, other.y, other.z]}')

@dataclass(init=False)
class Residue:
    """Class for storing Residues (defined as 3 atoms) belonging to a template
    and associated residue level information"""

    _atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]
    _vec: Vec3
    _indices: Tuple[int, int]

    def __init__(self, atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]):
        self._vec, self._indices = self.calc_residue_orientation(atoms)
        self._atoms = atoms

    @classmethod
    def calc_residue_orientation(cls, atoms: Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]) -> Tuple[Vec3, Tuple[int, int]]:
        # dictionary in which the vectors from start to finish are defined for each aminoacid type
        # the orientation vector is calculated differently for different aminoacid types
        vector_atom_type_dict = {'GLY': ('C','O'),
                            'ALA': ('CA','CB'),
                            'VAL': ('CA','CB'),
                            'LEU': ('CA','CB'),
                            'ILE': ('CA','CB'),
                            'MET': ('CG','SD'),
                            'PHE': ('CZ','mid'),
                            'TYR': ('CZ','OH'),
                            'TRP': ('CZ2','NE1'),
                            'CYS': ('CB','SG'),
                            'PRO': ('C','O'),
                            'SER': ('CB','OG'),
                            'THR': ('CB','OG1'),
                            'ASN': ('CG','OD1'),
                            'GLN': ('CD','OE1'),
                            'LYS': ('CE','NZ'),
                            'ARG': ('CZ','mid'),
                            'HIS': ('CG','ND1'),
                            'ASP': ('CG','mid'),
                            'GLU': ('CD','mid'),
                            'PTM': ('CA','CB'),
                            'ANY': ('C','O')}
        try:
            vectup = vector_atom_type_dict[atoms[0].residue_names[0]]
        except KeyError as exc:
                raise KeyError(f'Residue orientation is not defined for the residue type {atoms[0].residue_names[0]}') from exc

        # In residues with two identical atoms, the vector is calculated from the middle atom to the mid point between the identical pair
        if vectup[1] == 'mid':
            try:
                middle_index, middle_atom = next((index, atom) for index, atom in enumerate(atoms) if atom.atom_names[0] == vectup[0])
                side1, side2 = [atom for atom in atoms if atom != middle_atom]
                midpoint = (Vec3.from_xyz(side1) + Vec3.from_xyz(side2)) / 2
                return midpoint - Vec3.from_xyz(middle_atom), (middle_index, 9)
            except StopIteration:
                raise ValueError(f"Failed to find middle atom for amino-acid {atoms[0].residue_names[0]!r}") from None

        else: # from first atom to second atom
            try:
                first_atom_index, first_atom = next((index, atom) for index, atom in enumerate(atoms) if atom.atom_names[0] == vectup[0])
            except StopIteration:
                raise ValueError(f'Failed to find first atom for amino-acid {atoms[0].residue_names[0]!r}') from None
            try:
                second_atom_index, second_atom = next((index, atom) for index, atom in enumerate(atoms) if atom.atom_names[0] == vectup[1])
            except StopIteration:
                raise ValueError(f'Failed to find second atom for amino-acid {atoms[0].residue_names[0]!r}') from None
            return Vec3.from_xyz(second_atom) - Vec3.from_xyz(first_atom), (first_atom_index, second_atom_index)

    @property
    def atoms(self) -> Tuple[pyjess.TemplateAtom, pyjess.TemplateAtom, pyjess.TemplateAtom]:
        return self._atoms

    @property
    def residue_name(self) -> str:
        """`str`: Get the amino-acid type as three letter code from the first atom
        """
        return self.atoms[0].residue_names[0]

    @property
    def allowed_residues(self) -> str:
        """`str`: Get the allowed residue types as string of single letter codes
        """
        convert_to_single = {'ALA': 'A',
                            'CYS': 'C',
                            'ASP': 'D',
                            'GLU': 'E',
                            'PHE': 'F',
                            'GLY': 'G',
                            'HIS': 'H',
                            'ILE': 'I',
                            'LYS': 'K',
                            'LEU': 'L',
                            'MET': 'M',
                            'ASN': 'N',
                            'PRO': 'P',
                            'GLN': 'Q',
                            'ARG': 'R',
                            'SER': 'S',
                            'THR': 'T',
                            'VAL': 'V',
                            'TRP': 'W',
                            'TYR': 'Y',
                            'XXX': 'X'}
        return ''.join(set(convert_to_single[i] for i in self.atoms[0].residue_names))

    @property
    def match_mode(self) -> int:
        """`int`: Get the residue specificity (interger) from the first atom. <100 is specific
        """
        return self.atoms[0].match_mode

    @property
    def backbone(self) -> bool:
        """`bool`: Boolean, An atom is a backbone atom if it may match ANY or XXX in residue_names
        """
        return ('ANY' in self.atoms[0].residue_names or 'XXX' in self.atoms[0].residue_names)

    @property
    def residue_number(self) -> int:
        """`int`: Get the pdb residue number from the first atom
        """
        return self.atoms[0].residue_number

    @property
    def chain_id(self) -> str:
        """`str`: Get the pdb chain_id id from the first atom
        """
        return self.atoms[0].chain_id

    @property
    def orientation_vector(self) -> Vec3:
        """`Vec3`: Calculate the residue orientation vector according to the residue type
        """
        return self._vec

    @property
    def orientation_vector_indices(self) -> Tuple[int, int]:
        """`tuple`: Return the indices of the atoms 
        between which the orientation vector was calculated according to the residue type
        """
        return self._indices

@dataclass(frozen=True)
class Cluster:
    id: int
    member: int
    size: int

    def __post_init__(self):
        if self.member > self.size:
            raise ValueError("Cluster member cannot be greater than cluster size")

class AnnotatedTemplate(pyjess.Template):
    """Class for storing Templates and associated information. Inherits from pyjess.Template"""

    residues: Iterable[Residue] = ()

    _CATH_MAPPING: ClassVar[Dict[str, List[str]]]
    _EC_MAPPING: ClassVar[Dict[str, str]]
    _PDB_SIFTS: ClassVar[Dict[str, Dict[str, List[str]]]]

    pdb_id: Optional[str] = field(default=None)
    template_id_string : Optional[str] = field(default=None)
    mcsa_id: Optional[int] = field(default=None)
    cluster: Optional[Cluster] = field(default=None)
    uniprot_id: Optional[str] = field(default=None)
    organism: Optional[str] = field(default=None)
    organism_id: Optional[int] = field(default=None)
    resolution: Optional[float] = field(default=None)
    experimental_method: Optional[str] = field(default=None)
    enzyme_discription: Optional[str] = field(default=None)
    represented_sites: Optional[int] = field(default=None)
    ec: List[str] = field(default_factory=list)
    cath: List[str] = field(default_factory=list)
    
    def __init__(self,
            atoms: Sequence[pyjess.TemplateAtom],
            id: str,
            *,
            pdb_id: Optional[str] = None,
            template_id_string: Optional[str] = None,
            mcsa_id: Optional[int] = None,
            cluster: Optional[Cluster] = None,
            uniprot_id: Optional[str] = None,
            organism: Optional[str] = None,
            organism_id: Optional[str] = None, # TODO this is to avoid a bug when I get a list of organism_id's like 5,6
            resolution: Optional[float] = None,
            experimental_method: Optional[str] = None,
            enzyme_discription: Optional[str] = None,
            represented_sites: Optional[int] = None,
            ec: Iterable[str] = (),
            cath: Iterable[str] = ()
        ):
            super().__init__(atoms, id)

            residues = []
            for atom_triplet in chunks(atoms, 3): # yield chunks of 3 atom lines each
                if len(atom_triplet) != 3:
                    raise ValueError(f'Failed to construct residues. Got only {len(atom_triplet)} ATOM lines')
                # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
                unique_residues = { (atom.residue_names[0], atom.chain_id, atom.residue_number) for atom in atom_triplet }
                if len(unique_residues) != 1:
                    raise ValueError(f'Failed to construct residues. Atoms of different match_mode, chains, residue types or residue numbers {unique_residues} found in Atom triplet')
                residues.append(Residue(atom_triplet))

            super().__setattr__("residues", tuple(residues))
            super().__setattr__("pdb_id", pdb_id)
            super().__setattr__("mcsa_id", mcsa_id)
            super().__setattr__("template_id_string", template_id_string)
            super().__setattr__("cluster", cluster)
            super().__setattr__("uniprot_id", uniprot_id)
            super().__setattr__("organism", organism)
            super().__setattr__("organism_id", organism_id)
            super().__setattr__("resolution", resolution)
            super().__setattr__("experimental_method", experimental_method)
            super().__setattr__("enzyme_discription", enzyme_discription)
            super().__setattr__("represented_sites", represented_sites)

            super().__setattr__("ec", sorted({*ec, *self._add_ec_annotations()}))
            super().__setattr__("cath", sorted({*cath, *self._add_cath_annotations()}))

    @classmethod # reading from text
    def annotated_loads(cls, text: str, internal_template_id: str, warn: bool = True) -> AnnotatedTemplate:
        return cls.annotated_load(io.StringIO(text), internal_template_id, warn=warn)

    @classmethod # reading from TextIO
    def annotated_load(cls, file: TextIO, internal_template_id: str, warn: bool = True) -> AnnotatedTemplate:
        """Overloaded load to parse a pyjess.Template and its associated info into a AnnotatedTemplate object from TextIO
        """
        atoms = []
        metadata: dict[str, object] = {
            "ec": list(),
            "cath": list()
        }

        _PARSERS = {
            'ID' : cls._parse_template_id_string,
            'PDB_ID': cls._parse_pdb_id,
            'UNIPROT_ID': cls._parse_uniprot_id,
            'MCSA_ID': cls._parse_mcsa_id,
            'CLUSTER': cls._parse_cluster,
            'ORGANISM_NAME': cls._parse_organism_name,
            'ORGANISM_ID': cls._parse_organism_id,
            'RESOLUTION': cls._parse_resolution,
            'EXPERIMENTAL_METHOD': cls._parse_experimental_method,
            'EC': cls._parse_ec,
            'CATH': cls._parse_cath,
            'ENZYME': cls._parse_enzyme_discription,
            'REPRESENTING': cls._parse_represented_sites
        }
                
        # Note that currently mutliple templates per file are not supported by this parser function.
        # Note this parser does not check for the existance of a header like REMARK TEMPLATE
        atom_lines = []
        for line in filter(str.strip, file):
            tokens = line.split()
            if tokens[0] == 'REMARK':
                if warn and len(tokens) == 1:
                    warnings.warn('Expected some annotation information after the REMARK flag. Line Skipped.')
                    continue
                parser = _PARSERS.get(tokens[1])
                if parser is not None:
                    if len(tokens) < 3:
                        raise IndexError(f'Expected some annotation after the REMARK {tokens[1]} flag')
                    elif tokens[2].upper() in ['NONE', '?', 'NAN', 'NA']:
                        continue
                    parser(tokens, metadata, warn=warn)
            elif tokens[0] == 'ATOM':
                if line in atom_lines:
                    raise ValueError('Duplicate Atom lines passed!')
                atom_lines.append(line)
                atoms.append(pyjess.TemplateAtom.loads(line))
            elif tokens[0] == 'HETATM':
                raise ValueError('Supplied template with HETATM record. HETATMs cannot be searched by Jess')
            else:
                continue

        return AnnotatedTemplate(
            atoms = atoms,
            id = internal_template_id,
            **metadata, # type: ignore # unpack everything parsed into metadata
        )

    def dumps(self) -> str:
        buffer = io.StringIO()
        self.dump(buffer)
        return buffer.getvalue() # returns entire content temporary file object as a string

    def dump(self, file: IO[str]):
        file.write('REMARK TEMPLATE\n')
        if self.template_id_string:
            file.write(f'REMARK ID {self.template_id_string}\n')
        if self.pdb_id:
            file.write(f"REMARK PDB_ID {self.pdb_id}\n")
        if self.mcsa_id:
            file.write(f"REMARK MCSA_ID {self.mcsa_id}\n")
        if self.uniprot_id:
            file.write(f"REMARK UNIPROT_ID {self.uniprot_id}\n")
        if self.cluster: 
            file.write(f"REMARK CLUSTER {'_'.join([str(self.cluster.id), str(self.cluster.member), str(self.cluster.size)])}\n")
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
            file.write(f"REMARK REPRESENTING {self.represented_sites} CATALYTIC SITES\n")
        if self.experimental_method:
            file.write(f"REMARK EXPERIMENTAL_METHOD {self.experimental_method}\n")
        if self.ec:
            file.write(f"REMARK EC {','.join(self.ec)}\n")
        if self.cath:
            file.write(f"REMARK CATH {','.join(self.cath)}")
        file.write(f'REMARK SIZE {self.effective_size}\n')
        file.write(f'REMARK DIMENSION {self.dimension}\n')
        file.write(f'REMARK MULTIMERIC {self.multimeric}\n')
        file.write(f'REMARK RELATIVE_ORDER {self.relative_order}\n')

        for residue in self.residues:
            file.write(f'REMARK ORIENTATION_VECTOR OF RESIDUE {residue.residue_number}: between atom {residue.orientation_vector_indices[0]} and {residue.orientation_vector_indices[1]} {residue.orientation_vector.x:.3f} {residue.orientation_vector.y:.3f} {residue.orientation_vector.z:.3f}\n')

        for residue in self.residues:
            for atom in residue.atoms:
                raise NotImplementedError('Still need to write a dump method for templates in PyJess')
                # atom.dump(file)

        file.write('END\n')
    
    @cached_property
    def effective_size(self) -> int:
        """`int`: The number of unique residues in the template, excluding backbone residues and unspecific residues.
        """
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
            if residue.match_mode < 100 and not residue.backbone: # type specific and not Backbone
                effective_size += 1
        return effective_size

    @cached_property
    def multimeric(self) -> bool: # if the template is split across multiple protein chains
        """`bool`: Boolean if the template residues stem from multiple protein chains
        """
        return not all(res.chain_id == list(self.residues)[0].chain_id for res in self.residues)

    @cached_property
    def relative_order(self) -> List[int]: # list with length of deduplicated template dimension
        """`list` of `int`: Relative order of residues in the template sorted by the pdb residue number. This only works for non-multimeric templates
        """
        if self.multimeric:
            return [0]
        else:
            # Now extract relative template order
            return ranked_argsort([res.residue_number for res in self.residues])
    
    def _add_cath_annotations(self) -> List[str]:
        """`list` of `str`: Pull CATH Ids associated with that template from SIFTS and from the M-CSA"""
        cath_list = []
        if self.mcsa_id:
            cath_list.extend(self._CATH_MAPPING[str(self.mcsa_id)]) # this is a bit inaccurate possibly ... shouldnt we pull via pdb_id....
        if self.pdb_id:
            pdbchains = set()
            for res in self.residues:
                pdbchains.add("{}{}".format(self.pdb_id, res.chain_id))

            # Iterating of all pdbchains which were part of the AnnotatedTemplate
            for pdbchain in pdbchains:
                # also include CATH annotations from PDB-SIFTS
                subdict = self._PDB_SIFTS.get(pdbchain)
                if subdict is not None:
                    sifts_caths = subdict.get('cath').copy() # type: ignore
                    for cath in sifts_caths:
                        if cath != '?': # type: ignore
                            cath_list.append(cath) # type: ignore

        return cath_list

    def _add_ec_annotations(self) -> List[str]:
        """`list` of `str`: Pull EC Annotations associated with that template from SIFTS and from the M-CSA"""
        ec_list = []
        if self.mcsa_id is not None:
            ec_list.append(self._EC_MAPPING[str(self.mcsa_id)]) # this is a bit inaccurate possibly ... shouldnt we pull via pdb_id....
        if self.pdb_id is not None:
            pdbchains = set()
            for res in self.residues:
                pdbchains.add("{}{}".format(self.pdb_id, res.chain_id))

            # Iterating of all pdbchains which were part of the AnnotatedTemplate
            for pdbchain in pdbchains:
                # also include EC annotations from PDB-SIFTS
                subdict = self._PDB_SIFTS.get(pdbchain)
                if subdict is not None:
                    sifts_ecs = subdict.get('ec').copy() # type: ignore
                    for ec in sifts_ecs:
                        if ec != '?': # type: ignore
                            ec_list.append(ec) # type: ignore

        return ec_list

    @classmethod
    def _parse_pdb_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        if len(tokens[2]) != 4:
            raise ValueError(f'Found {tokens[2]} which has more than the expected 4 characters of a PDB_ID.')
        metadata['pdb_id'] = tokens[2].lower()

    @classmethod
    def _parse_template_id_string(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['template_id_string'] = tokens[2]

    @classmethod
    def _parse_uniprot_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', tokens[2])
        if match:
            metadata['uniprot_id'] = match.group()
        else:
            raise ValueError(f'Did not find a valid UniProt ID, found {tokens[2]}')

    @classmethod
    def _parse_mcsa_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['mcsa_id'] = int(tokens[2])
        except ValueError as exc:
            raise ValueError(f'Did not find a M-CSA ID, found {tokens[2]}') from exc

    @classmethod
    def _parse_cluster(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            cluster = Cluster(*list(map(int, tokens[2].split('_'))))
            metadata['cluster'] = cluster
        except ValueError as exc:
            raise ValueError(f'Did not find a Cluster specification in the form <id>_<member>_<size>, found {tokens[2]}') from exc

    @classmethod
    def _parse_organism_name(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['organism'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_organism_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['organism_id'] = tokens[2]

    @classmethod
    def _parse_resolution(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['resolution'] = float(tokens[2])
        except ValueError as exc:
            raise ValueError(f'Ill-formatted pdb resolution: {tokens[2]}') from exc

    @classmethod
    def _parse_experimental_method(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['experimental_method'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_ec(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        matches = [match.group() for match in re.finditer(r'[1-7](\.(\-|\d{1,})){3}', tokens[2])]
        non_cat_matches = [match.group() for match in re.finditer(r'[1-7](\.(\-|\d{1,}|n\d{1,})){3}', tokens[2])]
        if matches:
            for ec in matches:
                if ec not in metadata['ec']: # type: ignore
                    metadata['ec'].append(ec)  # type: ignore
        elif non_cat_matches:
            for ec in non_cat_matches:
                if ec not in metadata['ec']: # type: ignore
                    metadata['ec'].append(ec)  # type: ignore
            if warn:
                warnings.warn(f'Rare EC number(s) {[ec for ec in non_cat_matches]} presumed to be noncatalytic detected!')
        else:
            raise ValueError(f'Did not find a valid EC number, found {tokens[2]}')

    @classmethod
    def _parse_cath(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        matches = [match.group() for match in re.finditer(r'[1-46](\.(\-|\d{1,})){3}', tokens[2])]
        if matches:
            for cath in matches:
                if cath not in metadata['cath']: # type: ignore
                    metadata['cath'].append(cath)  # type: ignore
        else:
            raise ValueError(f'Did not find a valid CATH number, found {tokens[2]}')

    @classmethod
    def _parse_enzyme_discription(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['enzyme_discription'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_represented_sites(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['represented_sites'] = int(tokens[2])
        except ValueError as exc:
            raise ValueError(f'Ill-formatted number of represented sites: {tokens[2]}') from exc

# Populate the mapping of MCSA IDs to CATH numbers so that it can be accessed
# by individual templates in the `Template.cath` property.
# Source: M-CSA which provides cath annotations for either residue homologs or for m-csa entries
with files(__package__).joinpath('data', 'MCSA_CATH_mapping.json').open() as f: #type: ignore
    AnnotatedTemplate._CATH_MAPPING = json.load(f)

with files(__package__).joinpath('data', 'MCSA_EC_mapping.json').open() as f: #type: ignore
    AnnotatedTemplate._EC_MAPPING = json.load(f)

# global MCSA_interpro_dict
# # dictonariy mapping M-CSA entries to Interpro Identifiers
# # Interpro Acceccesions at the Domain, Family and Superfamily level
# # are searched for the reference sequences of each M-CSA entry.
# # Note that an M-CSA entry may have multiple reference sequences

# Source: CATH, EC and InterPro from PDB-SIFTS through mapping to the pdbchain
with files(__package__).joinpath('data', 'pdb_sifts.json').open() as f: #type: ignore
    AnnotatedTemplate._PDB_SIFTS = json.load(f)

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
        raise FileNotFoundError(f'No template files with the {extension} extension found in the {directory_path.resolve()} directory')

def load_templates(template_dir: Path = Path(str(files(__package__).joinpath('jess_templates_20230210'))), warn: bool = True, verbose: bool = False, cpus: int = 0) -> Iterator[AnnotatedTemplate]:
    """Load templates from a given directory, recursively.
    """
    if not template_dir.exists():
        raise FileNotFoundError(template_dir)
    elif not template_dir.is_dir():
        raise NotADirectoryError(template_dir)
    
    template_paths = _get_paths_by_extension(template_dir, '.pdb')

    if verbose:
        print(f'Loading Template files from {str(template_dir.resolve())}')

    
    def _load_and_annotate(template_path, warn, template_index):
        try:
            with template_path.open() as f:
                return AnnotatedTemplate.annotated_load(file=f,
                                                        warn=warn,
                                                        internal_template_id=str(template_index)
                                                        )
        except ValueError as exc:
            raise ValueError(f'Passed Template file {template_path.resolve()} contained ATOM lines which are not in Jess Template format.') from exc

    if cpus == 0:
        cpus = os.cpu_count() or 1
    elif cpus < 0:
        cpus = max(1, (os.cpu_count() or 1) - cpus)

    pool = DummyPool() if cpus == 1 else ThreadPool(cpus)
    with pool:
            args = [(path, warn, idx) for idx, path in enumerate(template_paths)]
            try:
                template_iterator = pool.starmap(_load_and_annotate, args)
                for loaded_template in template_iterator:
                    yield loaded_template
            except ValueError as exc:
                raise exc

def check_template(template: AnnotatedTemplate, warn: bool = True) -> bool:
    # TODO improve this and write tests
    if warn:
        checker = True

        # Raise warnings if some properties could not be annotated!
        if not template.ec:
            checker = False
            warnings.warn(f'Could not find EC number annotations for the template file {template.id}')

        if not template.cath:
            checker = False
            warnings.warn(f'Could not find CATH annotations for the following template file {template.id}')

        if template.pdb_id:
            # check overlap between sifts mapping and CATH, EC annotations
            # Source: pdb to sifts mapping which maps CATH to pdb chain IDs and UniProt IDs, sifts also provides UniProt to EC mapping
            # Since AnnotatedTemplates may contain residues from multiple chains
            pdbchains = set()
            for res in template.residues:
                pdbchains.add(template.pdb_id + res.chain_id)

            # Iterating of all pdbchains which were part of the AnnotatedTemplate
            for pdbchain in pdbchains:
                # also include EC annotations from PDB-SIFTS
                subdict = template._PDB_SIFTS.get(pdbchain, None)
                if subdict:
                    # if template.cath and subdict['cath']:
                    #     if not set(template.cath) & set(subdict['cath']):
                    #         warnings.warn(f'Did not find an intersection of CATH domains as annotated by the M-CSA ID {template.mcsa_id} with {template.cath} versus PDF-SIFTS {template._PDB_SIFTS[pdbchain]['cath']} for template {filepath} from PDB ID {template.pdb_id}')
                    sifts_uniprot = subdict['uniprot_id']
                    if template.uniprot_id != sifts_uniprot:
                        checker = False
                        warnings.warn(f'Different UniProt Accessions {template.uniprot_id} and {sifts_uniprot} found in template {template.id} from PDB ID {template.pdb_id}')

        return checker
    else:
        return True
