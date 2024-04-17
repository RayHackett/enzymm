from __future__ import annotations

import glob
import warnings
import re
import json
from importlib.resources import files
from pathlib import Path
from typing import Sequence, List, Tuple, Set, Dict, Optional, Union, TextIO, Iterator, ClassVar, IO
from functools import cached_property
from dataclasses import dataclass, field
import io
import tempfile

import pyjess # type: ignore

from .utils import chunks, ranked_argsort

__all__ = [
    "Vec3",
    "Residue",
    "Atom",
    "Template",
    "load_templates",
]

# Find all template files which end in .pdb
def get_template_paths(directory_path: Path, extension: str ='.pdb') -> List[Path]:
    pattern = f"{directory_path}/**/*{extension}"
    file_paths = glob.glob(pattern, recursive=True)
    files = []
    for file in file_paths:
        files.append(Path(file))

    if files:
        return files
    else:
        raise FileNotFoundError(f'No template files with the {extension} extension found in the {directory_path.resolve()} directory')

@dataclass
class Vec3:
    '''Class for storing 3D vectors in XYZ'''
    x: float
    y: float
    z: float

@dataclass
class Residue:
    """Class for storing Residues (defined as 3 atoms) belonging to a Template and associated residue level information"""
    atoms: Tuple[Atom, Atom, Atom]

    def __post_init__(self) -> None:
        if not self.atoms:
            raise ValueError("Cannot create a residue with no atoms")

    @classmethod
    def calc_residue_orientation(cls, atoms: Tuple[Atom, Atom, Atom]) -> Vec3:
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
        
        vectup = vector_atom_type_dict[atoms[0].residue_name]
        # In residues with two identical atoms, the vector is calculated between the middle atom and the mid point between the identical pair
        if vectup[1] == 'mid':
            try:
                middle = next(Atom for Atom in atoms if Atom.name == vectup[0])
                side1, side2 = [atom for atom in atoms if atom != middle]
                midpoint = [(side1.x + side2.x) / 2, (side1.y + side2.y) / 2, (side1.z + side2.z) / 2]
                return Vec3(middle.x - midpoint[0], middle.y - midpoint[1], middle.z - midpoint[2])
            except StopIteration:
                raise ValueError(f"Failed to find middle atom for amino-acid {atoms[0].residue_name!r}") from None

        else: # from first atom to second atom
            try:
                first_atom = next(Atom for Atom in atoms if Atom.name == vectup[0])
            except StopIteration:
                raise ValueError(f'Failed to find first atom for amino-acid {atoms[0].residue_name!r}') from None
            try:
                second_atom = next(Atom for Atom in atoms if Atom.name == vectup[1])
            except StopIteration:
                raise ValueError(f'Failed to find second atom for amino-acid {atoms[0].residue_name!r}') from None
                
            return first_atom - second_atom

    @property
    def residue_name(self) -> str:
        """`str`: Get the amino-acid type as three letter code from the first atom
        """
        return self.atoms[0].residue_name

    @property
    def allowed_residues(self) -> str:
        """`str`: Get the allowed residue types as string of single letter codes from the first atom
        """
        return self.atoms[0].allowed_residues

    @property
    def specificity(self) -> int:
        """`int`: Get the residue specificity (interger) from the first atom. <100 is specific
        """
        return self.atoms[0].specificity

    @property
    def backbone(self) -> bool:
        """`int`: Get boolean indicated if the residue is composed of backbone atoms from the first atom
        """
        return self.atoms[0].backbone

    @property
    def residue_number(self) -> int:
        """`int`: Get the pdb residue number from the first atom
        """
        return int(self.atoms[0].residue_number)

    @property
    def chain_id(self) -> str:
        """`int`: Get the pdb chain_id id from the first atom
        """
        return self.atoms[0].chain_id

    @property
    def orientation_vector(self) -> Vec3:
        """`int`: Calculate the residue orientation vector according to the residue type
        """
        return self.calc_residue_orientation(self.atoms)

@dataclass
class Atom:
    '''Class for storing Atom information from a Template file'''
    specificity: int
    name: str
    residue_name: str
    backbone: bool
    chain_id: str
    residue_number: int
    x: float
    y: float
    z: float
    allowed_residues: str
    altloc: Optional[str] = field(default=None)
    flexibility: Optional[float] = field(default=None)

    @classmethod
    def loads(cls, line: str) -> Atom:
        ## this uses a modified PDB format
        # first four characters are 'ATOM' [0:4]
        specificity = int(line[7:11]) # inplace of atom number
        name = str(line[12:16].strip()) # atom_name
        altloc = str(line[16].strip()) # no idea what this is used for in a template, optional
        residue_name = str(line[17:20].strip())
        backbone = line[17:20].strip() == 'ANY'
        chain_id = str(line[20:22].strip()) # allows for 2 letter chain IDs
        residue_number = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        allowed_residues = str(line[55:60].strip()) # inplace of Occupancy
        flexibility = float(line[61:66]) # inplace of Temperature Factor, optional
        # segment identifier: empty
        # element: empty
        # charge: empty
        return cls(specificity=specificity, name=name, altloc=altloc, residue_name=residue_name, backbone=backbone, chain_id=chain_id, residue_number=residue_number, x=x, y=y, z=z, allowed_residues=allowed_residues, flexibility=flexibility)

    def dumps(self):
        buffer = io.StringIO()
        self.dump(buffer)
        return buffer.getvalue() # returns entire content temporary file object as a string

    def dump(self, file: IO[str]):
        if self.flexibility:
            pdb_line = f"ATOM  {self.specificity:>5} {self.name:^4}{self.altloc if self.altloc is not None else '':<1}{self.residue_name:<3}{self.chain_id:>2}{self.residue_number:>4}    {self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f} {self.allowed_residues:<5}{self.flexibility:>5.2f} \n"
        else:
            pdb_line = f"ATOM  {self.specificity:>5} {self.name:^4}{self.altloc if self.altloc is not None else '':<1}{self.residue_name:<3}{self.chain_id:>2}{self.residue_number:>4}    {self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f} {self.allowed_residues:<5}{'':<5} \n"
        file.write(pdb_line)

    def __sub__(self, other):
        if not isinstance(other, Atom):
            raise TypeError(f'Expected an instance of "Atom" but got {type(other)}')
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

@dataclass
class Template:
    """Class for storing Templates and associated information"""
    residues: List[Residue]
    pdb_id: str
    mcsa_id: int
    cluster: str
    uniprot_id: Optional[str] = field(default=None)
    organism: Optional[str] = field(default=None)
    organism_id: Optional[int] = field(default=None)
    resolution: Optional[float] = field(default=None)
    experimental_method: Optional[str] = field(default=None)
    enzyme_discription: Optional[str] = field(default=None)
    represented_sites: Optional[int] = field(default=None)
    ec: Optional[Set[str]] = field(default=None)

    _CATH_MAPPING: ClassVar[Dict[str, List[str]]]
    _EC_MAPPING: ClassVar[Dict[str, List[str]]]
    _PDB_SIFTS: ClassVar[Dict[str, Dict[str, List[str]]]]

    @classmethod # reading from text
    def loads(cls, text: str, warn: bool = True) -> Template:
        return cls.load(io.StringIO(text), warn=warn)

    @classmethod # reading from TextIO
    def load(cls, file: TextIO, warn: bool = True) -> Template:
        """Function to parse a template and its associated info into a Template object from TextIO
        """
        atom_lines = []
        residues = []
        metadata: dict[str, object] = {
            "ec": set(),
        }

        _PARSERS = {
            'PDB_ID': cls._parse_pdb_id,
            'UNIPROT_ID': cls._parse_uniprot_id,
            'MCSA_ID': cls._parse_mcsa_id,
            'CLUSTER': cls._parse_cluster,
            'ORGANISM_NAME': cls._parse_organism_name,
            'ORGANISM_ID': cls._parse_organism_id,
            'RESOLUTION': cls._parse_resolution,
            'EXPERIMENTAL_METHOD': cls._parse_experimental_method,
            'EC': cls._parse_ec,
            'ENZYME': cls._parse_enzyme_discription,
            'REPRESENTING': cls._parse_represented_sites
        }
                
        for line in filter(str.strip, file):
            tokens = line.split()
            if tokens[0] == 'REMARK':
                parser = _PARSERS.get(tokens[1])
                if parser is not None:
                    parser(tokens, metadata, warn=warn)
            elif tokens[0] == 'ATOM':
                atom_lines.append(line)
            else:
                continue

        for atom_triplet in chunks(atom_lines, 3): # yield chunks of 3 atom lines each
            atoms: Tuple[Atom, Atom, Atom] = tuple(Atom.loads(line) for line in atom_triplet)  # type: ignore
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = { (atom.residue_name, atom.chain_id, atom.residue_number) for atom in atoms }
            if len(unique_residues) != 1:
                raise ValueError('Multiple residues found in atom lines')
            residues.append(Residue(atoms))

        # also include EC annotations from the M-CSA mapping
        metadata['ec'].add(cls._EC_MAPPING[str(metadata['mcsa_id'])])  # type: ignore

        # Since Templates may contain residues from multiple chains
        pdbchains = set()
        for res in residues:
            pdbchains.add("{}{}".format(metadata['pdb_id'], res.chain_id))

        # Iterating of all pdbchains which were part of the Template
        for pdbchain in pdbchains:
            # also include EC annotations from PDB-SIFTS
            subdict = cls._PDB_SIFTS.get(pdbchain)
            if subdict is not None: # TODO this is an ugly workaround: Technically this is a bug either with missing annotations in SIFTS or due to wierd chain name conventions in .cif files
                metadata['ec'].update(subdict.get('ec').copy())  # type: ignore

        return Template(
            residues = residues,
            **metadata, # type: ignore # unpack everything we parsed into metadata
        )

    def dumps(self):
        buffer = io.StringIO()
        self.dump(buffer)
        return buffer.getvalue() # returns entire content temporary file object as a string

    def dump(self, file: IO[str]):
        file.write(f'REMARK TEMPLATE\n')
        file.write(f'REMARK ID {self.id}\n')
        file.write(f'REMARK PDB_ID {self.pdb_id}\n')
        file.write(f'REMARK MCSA_ID {self.mcsa_id}\n')
        file.write(f'REMARK UNIPROT_ID {self.uniprot_id}\n')
        file.write(f'REMARK CLUSTER {self.cluster}\n')
        file.write(f'REMARK ORGANISM_NAME {self.organism}\n')
        file.write(f'REMARK ORGANISM_ID {self.organism_id}\n')
        file.write(f'REMARK RESOLUTION {self.resolution}\n')
        file.write(f'REMARK EXPERIMENTAL_METHOD {self.experimental_method}\n')
        file.write(f'REMARK ENZYME {self.enzyme_discription}\n')
        file.write(f'REMARK REPRESENTING {self.represented_sites} CATALYTIC SITES\n')
        file.write(f'REMARK EXPERIMENTAL_METHOD {self.experimental_method}')
        file.write(f"REMARK EC {','.join(self.ec) if self.ec is not None else ''}\n")
        file.write(f"REMARK CATH {','.join(self.cath) if self.cath else ''}\n")
        file.write(f'REMARK SIZE {self.size}\n')
        file.write(f'REMARK TRUE_SIZE {self.true_size}\n')
        file.write(f'REMARK MULTIMERIC {self.multimeric}\n')
        file.write(f'REMARK RELATIVE_ORDER {self.relative_order}\n')

        for residue in self.residues:
            file.write(f'REMARK ORIENTATION_VECTOR OF RESIDUE {residue.residue_number}: {residue.orientation_vector.x:.3f} {residue.orientation_vector.y:.3f} {residue.orientation_vector.z:.3f}\n')

        for residue in self.residues:
            for atom in residue.atoms:
                atom.dump(file)

        file.write('END\n')


    def to_pyjess_template(self) -> pyjess.Template:
        return pyjess.Template.loads(self.dumps(), id=self.id)
        # with tempfile.NamedTemporaryFile(mode='w+') as f: # TODO is this bad? default is mode w+b
        #     self.dump(f) # writing to tempfile f
        #     f.flush() # makes sure the file is written to the disk and not in IO-buffer
        #     f.seek(0)
        #     return pyjess.Template.loads(f.read(), id=self.id)

    @property
    def id(self) -> str: # pyjess templates need a name, so mine do too
        """`str`: The name of the template composed of M-CSA ID, true size, pdb_id, cluster seperated by underscores
        """
        return f'{self.mcsa_id}_{self.true_size}_{self.pdb_id}_{self.cluster}'

    @property
    def true_size(self) -> int: # Number of residues in the template as counted by triplets of atoms
        """`int`: The number of residues in the template, as counted by triplets of atoms.
        """
        return len(self.residues)
    
    @cached_property
    def size(self) -> int:
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
            if residue.specificity < 100 and residue.backbone == False: # type specific and not Backbone
                effective_size += 1
        return effective_size

    @cached_property
    def multimeric(self) -> bool: # if the template is split across multiple protein chains
        """`bool`: Boolean if the template residues stem from multiple protein chains
        """
        return all(res.chain_id == self.residues[0].chain_id for res in self.residues)

    @cached_property
    def relative_order(self) -> List[int]: # list with length of deduplicated true_size
        """`list` of `int`: Relative order of residues in the Template sorted by the pdb residue number. This only works for non-multimeric Templates
        """
        if self.multimeric:
            return [0]
        else:
            # Now extract relative template order
            return ranked_argsort([res.residue_number for res in self.residues])
    
    @property
    def cath(self) -> Set[str]:
        """`set` of `str`: Set of CATH Ids associated with that Template"""
        return set(self._CATH_MAPPING[str(self.mcsa_id)].copy())

    @classmethod
    def _parse_pdb_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        if not tokens[2]:
            raise ValueError(f'Did not find a PDB ID, found {tokens[2]}')
        else:
            metadata['pdb_id'] = tokens[2].lower()

    @classmethod
    def _parse_uniprot_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', tokens[2])
        if match:
            metadata['uniprot_id'] = match.group()
        else:
            if warn:
                warnings.warn(f'Did not find a valid UniProt ID, found {tokens[2]}')

    @classmethod
    def _parse_mcsa_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['mcsa_id'] = int(tokens[2])
        except ValueError as exc:
                raise ValueError(f'Did not find a M-CSA ID, found {tokens[2]}') from exc

    @classmethod
    def _parse_cluster(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        if not tokens[2]:
            raise ValueError(f'Did not find a Cluster ID, found {tokens[2]}')
        else:
            metadata['cluster'] = tokens[2]

    @classmethod
    def _parse_organism_name(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['organism'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_organism_id(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['organism_id'] = int(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Ill-formatted organism ID: {tokens[2]}')

    @classmethod
    def _parse_resolution(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['resolution'] = float(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Ill-formatted pdb resolution: {tokens[2]}')

    @classmethod
    def _parse_experimental_method(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['experimental_method'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_ec(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        matches = [match.group() for match in re.finditer(r'\d{1,2}(\.(\-|\d{1,2})){3}', tokens[2])]
        non_cat_matches = [match.group() for match in re.finditer(r'\d{1,2}(\.(\-|\d{1,2}|n\d{1,2})){3}', tokens[2])] 
        if matches:
            metadata['ec'].update(matches)  # type: ignore
        elif non_cat_matches:
            metadata['ec'].update(non_cat_matches)  # type: ignore
            if warn:
                warnings.warn(f'Rare EC number {non_cat_matches} presumed to be noncatalytic detected!')
        else:
            if warn:
                warnings.warn(f'Did not find a valid EC number, found {tokens[2]}')

    @classmethod
    def _parse_enzyme_discription(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        metadata['enzyme_discription'] = ' '.join(tokens[2:])

    @classmethod
    def _parse_represented_sites(cls, tokens: List[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['represented_sites'] = int(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Ill-formatted number of represented sites: {tokens[2]}')

# Populate the mapping of MCSA IDs to CATH numbers so that it can be accessed
# by individual templates in the `Template.cath` property.
# Source: M-CSA which provides cath annotations for either residue homologs or for m-csa entries
with files(__package__).joinpath('data', 'MCSA_CATH_mapping.json').open() as f:
    Template._CATH_MAPPING = json.load(f)

with files(__package__).joinpath('data', 'MCSA_EC_mapping.json').open() as f:
    Template._EC_MAPPING = json.load(f)

# global MCSA_interpro_dict
# # dictonariy mapping M-CSA entries to Interpro Identifiers
# # Interpro Acceccesions at the Domain, Family and Superfamily level
# # are searched for the reference sequences of each M-CSA entry.
# # Note that an M-CSA entry may have multiple reference sequences

# Source: CATH, EC and InterPro from PDB-SIFTS through mapping to the pdbchain
with files(__package__).joinpath('data', 'pdb_sifts.json').open() as f:
    Template._PDB_SIFTS = json.load(f)

# # add a list of cofactors associated with each EC number from Neeras List
# cofactors = set()
# if Template_EC in cofactor_dict:
#     cofactors.update(cofactor_dict[Template_EC])

def load_templates(template_dir: Path = Path(str(files(__package__).joinpath('jess_templates_20230210'))), warn: bool = True) -> Iterator[Tuple[Template, Path]]:
    """Load templates from a given directory, recursively.
    """
    if not template_dir.exists():
        raise FileNotFoundError(template_dir)
    elif not template_dir.is_dir():
        raise NotADirectoryError(template_dir)

    template_paths = get_template_paths(template_dir)
    for template_path in template_paths:
        with template_path.open() as f:
            try:
                # Yield the Template and its filepath as a tuple
                yield (Template.load(file=f, warn=warn), template_path)
            except ValueError as exc:
                raise ValueError(f"Failed to parse {template_path}!") # from exc

def group_templates_by_size(templates: List[Template]) -> Dict[int, List[Template]]:
    grouped_templates: Dict[int, List[Template]] = {}
    for template in templates:
        if template.size not in grouped_templates:
            grouped_templates[template.size] = []
        grouped_templates[template.size].append(template)
    return grouped_templates

def check_template(Template_tuple: Tuple[Template, Path], warn: bool = True) -> bool:
    if warn:
        Template, filepath = Template_tuple

        # Raise warnings if some properties could not be annotated!
        if not Template.ec:
            warnings.warn(f'Could not find EC number annotations for the template file {filepath}')

        if not Template.cath:
            warnings.warn(f'Could not find CATH annotations for the following template file {filepath}')

        # check overlap between sifts mapping and CATH, EC annotations
        # Source: pdb to sifts mapping which maps CATH to pdb chain IDs and UniProt IDs, sifts also provides UniProt to EC mapping
        # Since Templates may contain residues from multiple chains
        pdbchains = set()
        for res in Template.residues:
            pdbchains.add(Template.pdb_id + res.chain_id)

        # Iterating of all pdbchains which were part of the Template
        for pdbchain in pdbchains:
            # also include EC annotations from PDB-SIFTS
            subdict = Template._PDB_SIFTS.get(pdbchain, None)
            if subdict: # TODO this is an ugly workaround: Technically this is a bug either with missing annotations in SIFTS or due to wierd chain name conventions in .cif files
                # if Template.cath and subdict['cath']:
                #     if not set(Template.cath) & set(subdict['cath']):
                #         warnings.warn(f'Did not find an intersection of CATH domains as annotated by the M-CSA ID {Template.mcsa_id} with {Template.cath} versus PDF-SIFTS {Template._PDB_SIFTS[pdbchain]['cath']} for Template {filepath} from pdb {Template.pdb_id}')
                sifts_uniprot = subdict['uniprot_id']
                if Template.uniprot_id != sifts_uniprot:
                    warnings.warn(f'Different UniProt Accessions {Template.uniprot_id} and {sifts_uniprot} found in Template {filepath} from pdbid {Template.pdb_id}')

        return True # change this to false if you want to exclude templates based on any of their properties!
    else:
        return True

if __name__ == "__main__":
    # TODO this is here for debugging purposes - delete later
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        templates = list(load_templates(Path(str(files(__package__).joinpath('jess_templates_20230210', f'{template_res_num}_residues'))), warn=False))
        print(f'current template size folder: {template_res_num}')
        print(len(templates))
        #for template in templates:
        #    check_template(template)

