from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import OrderedDict
import json
from typing import Optional, List, Dict, Tuple

# NOTE I know these could be frozen dataclasses too
# but that makes it annoying to apply manually fixes post init
# which is unfortunately necessary sometimes


@dataclass(frozen=False)
class HomologousResidue:
    """
    Base class for a PDB residue identical or homologous to a template residue.

    Attributes:
        name : `str` Three letter amino-acid code
        number : `int` residue id in the PDB
        auth_number : `int` Author assigned residue id

    """

    name: str
    number: int
    auth_number: int
    # domain_name: Optional[str]
    # domain_cath_id: Optional[str]

    def __post_init__(self):
        object.__setattr__(self, "name", self.name.upper())


@dataclass(frozen=False)
class ReferenceCatalyticResidue(HomologousResidue):
    """
    Child Class of `HomologousResidue` for a reference PDB residue for a template residue.

    Attributes:
        function_location_abv : `str` Functional part of the residue. Empty string means side chain. Else "main", "main-C", "main-N" or "ptm"
        ptm : `str` Empty string or for post translationally modified residues, the 3-letter code of the ptm in the PDB
        roles : `list` of `str` of EMO codes describing functional roles
        roles_summary : `list` of `str` discriptions of functional roles

    """

    function_location_abv: Optional[str]
    ptm: Optional[str]
    roles: List[str]
    roles_summary: List[str]


@dataclass(frozen=False)
class NonReferenceCatalyticResidue(HomologousResidue):
    """
    Child Class of `HomologousResidue` for a non-reference PDB residue for a template residue.

    Attributes:
        reference : `tuple` of (`int` , `str` , `int`) : Tuple of mcsa_id, pdb_id, residue_index to find the reference

    """

    reference: Tuple[int, str, int]  # can link to with mcsa_id, pdb_id, residue index


@dataclass(frozen=False)
class HomologousPDB:
    """
    Class for a PDB Entry homologous to or identical to the structure from which a template was generated.

    Attributes:
        mcsa_id : `int` : M-CSA ID
        reference_pdbchain : `str` : PDB identifier of the M-CSA reference PDB structure
        is_reference : `bool` : If this structure is an M-CSA reference
        pdb_id : `str` : PDB identifier of the PDB structure
        chain_id : `str` : Chain identifier in the PDB
        assembly_chain_id : `str` : Author assigned assembly chain identifier
        assembly : `int` : Assembly number (usually but not always the biological assembly)
        residues : `dict` of [`int` , `HomologousResidue`] : Dictionary mapping residue indices to residues. If "is_reference" then values are `ReferenceCatalyticResidue` else they are `NonReferenceCatalyticResidue`

    """

    # NOTE
    # I purposefully left out UniProt Identifers because for some PDBchains are
    # not easily assigned to a UniProt Identifer. Sifts does a great job but
    # some chains are not in UniProt / are not mapped
    # even worse, some chains have multiple covalently attached UniProts
    # and are chimeric (sometimes even across different species)

    # because we iteratively add residues, the HomologousPDB dataclass
    # cannot be frozen

    mcsa_id: int
    reference_pdbchain: str
    is_reference: bool
    pdb_id: str
    chain_id: str
    assembly_chain_id: str
    assembly: int
    residues: Dict[int, HomologousResidue]  # can fully contain

    def __post_init__(self):
        object.__setattr__(self, "pdb_id", self.pdb_id.lower())

    def _state(self) -> Tuple:
        return (
            self.mcsa_id,
            self.reference_pdbchain,
            self.is_reference,
            self.pdb_id,
            self.chain_id,
            self.assembly_chain_id,
            self.assembly,
            tuple(self.residues),
        )

    def __eq__(self, other: object) -> bool:
        if isinstance(other, HomologousPDB):
            self_state = self._state()
            other_state = other._state()
            return self_state == other_state
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, HomologousPDB):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self), self._state()))


# @dataclass(frozen=True)
# class Entry:
#     uniprot: str
#     reference_pdb: HomologousPDB
#     number_of_residues: int


def load_mcsa_catalytic_residue_homologs_info(
    outdir: Path,
) -> Dict[str, Dict[str, HomologousPDB]]:

    file = outdir.joinpath("catalytic_residue_homologs_information.json")
    # file is probably in Path(outdir, "catalytic_residue_homologs_information.json")
    with open(file, "r") as f:
        data = json.load(f)

    out: OrderedDict = OrderedDict()
    for mcsa_id, pdb_entries in data.items():
        mcsa_id = int(mcsa_id)
        out[mcsa_id] = OrderedDict()
        for pdb_id, pdb_data in pdb_entries.items():
            residues: Dict[int, HomologousResidue] = {}
            # If it is a reference structure,
            # for multimetic catalytic sites, references may be mixed
            if pdb_data["is_reference"]:
                for index, residue in pdb_data["residues"].items():
                    index = int(index)
                    try:
                        residues[index] = ReferenceCatalyticResidue(
                            name=residue["code"],
                            number=residue["resid"],
                            auth_number=residue["auth_resid"],
                            function_location_abv=residue["function_location_abv"],
                            ptm=residue["ptm"],
                            roles=residue["roles"],
                            roles_summary=residue["roles_summary"],
                        )
                    except KeyError:
                        residues[index] = NonReferenceCatalyticResidue(
                            name=residue["code"],
                            number=residue["resid"],
                            auth_number=residue["auth_resid"],
                            reference=tuple(residue["reference"]),
                        )
            # If it is not a reference, all residues are NonReference
            else:
                for index, residue in pdb_data["residues"].items():
                    index = int(index)
                    residues[index] = NonReferenceCatalyticResidue(
                        name=residue["code"],
                        number=residue["resid"],
                        auth_number=residue["auth_resid"],
                        reference=tuple(residue["reference"]),
                    )
            out[mcsa_id][pdb_id] = HomologousPDB(
                mcsa_id=mcsa_id,
                reference_pdbchain=pdb_data["reference_pdbchain"],
                is_reference=pdb_data["is_reference"],
                pdb_id=pdb_data["pdb_id"],
                chain_id=pdb_data["chain_name"],
                assembly_chain_id=pdb_data["assembly_chain_name"],
                assembly=pdb_data["assembly"],
                residues=residues,
            )

    return out
