from __future__ import annotations

import glob
import warnings
import re
from importlib.resources import files
from pathlib import Path
from typing import Optional, Union, TextIO, Iterator
from functools import cached_property
from dataclasses import dataclass, field

__all__ = [
    "Vec3",
    "Residue",
    "Atom",
    "Template",
    "load_templates",
]

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Find all template files which end in .pdb
def get_template_files(directory_path, extension='.pdb'):
    pattern = f"{directory_path}/**/*{extension}"
    file_names = glob.glob(pattern, recursive=True)
    files = []
    for file in file_names:
        files.append(Path(file))

    if files:
        return files
    else:
        raise ValueError('No template files with the .pdb extension found in this directory')

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

    # TODO TypeError: object.__init__() takes exactly one argument (the instance to initialize)
    def __init__(self, atoms: Tuple[Atom, Atom, Atom]) -> None:
        if not atoms:
            raise ValueError("Cannot create a residue with no atoms")
        super().__init__(atoms)

    @classmethod
    def get_vec_template(atoms: Tuple[Atom, Atom, Atom]) -> Vec3:
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
        
        vectup = vector_atom_type_dict[atoms[0].aa_type]
        # In residues with two identical atoms, the vector is calculated between the middle atom and the mid point between the identical pair
        if vectup[1] == 'mid':
            try:
                middle = next(Atom for Atom in atoms if Atom.atom_name == vectup[0])
                side1, side2 = [atom for atom in atoms if atom != middle]
                midpoint = [(side1.x + side2.x) / 2, (side1.y + side2.y) / 2, (side1.z + side2.z) / 2]
                return Vec3(middle.x - midpoint[0], middle.y - midpoint[1], middle.z - midpoint[2])
            except StopIteration:
                raise ValueError(f"Failed to find middle atom for amino-acid {atoms[0].aa_type!r}") from None

        else: # from first atom to second atom
            try:
                first_atom = next(Atom for Atom in atoms if Atom.atom_name == vectup[0])
            except StopIteration:
                raise ValueError(f'Failed to find first atom for amino-acid {atoms[0].aa_type!r}') from None
            try:
                second_atom = next(Atom for Atom in atoms if Atom.atom_name == vectup[1])
            except StopIteration:
                raise ValueError(f'Failed to find second atom for amino-acid {atoms[0].aa_type!r}') from None
                
            return first_atom - second_atom

    @property
    def aa_type(self) -> str:
        return self.atoms[0].aa_type

    @property
    def allowed_residues(self) -> str:
        return self.atoms[0].allowed_residues

    @property
    def specificity(self) -> int:
        return self.atoms[0].specificity

    @property
    def backbone(self) -> bool:
        return self.atoms[0].backbone

    @property
    def pdb_residue(self) -> int:
        return self.atoms[0].pdb_residue

    @property
    def chain(self) -> str:
        return self.atoms[0].chain

    @property
    def orientation_vector(self, atoms) -> Vec3:
        return self.get_vec_template(atoms)

@dataclass
class Atom:
    ''' Class for storing Atom information from a Template file'''
    x: float
    y: float
    z: float
    atom_name: str
    aa_type: str
    allowed_residues: str
    specificity: int
    backbone: bool
    pdb_residue: str
    chain: str

    @classmethod
    def loads(cls, line: str) -> Atom:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom_name = str(line[12:16].strip())
        aa_type = str(line[17:20].strip())
        allowed_residues = str(line[55:61].strip())
        specificity = int(line[8:11])
        backbone = line[17:20].strip() == 'ANY'
        pdb_residue = int(line[22:26])
        chain = str(line[20:22].strip())
        return cls(x=x, y=y, z=z, atom_name=atom_name, aa_type=aa_type, allowed_residues=allowed_residues, specificity=specificity, backbone=backbone, pdb_residue=pdb_residue, chain=chain)

    def __sub__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented #TODO is this a type error then?
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

@dataclass
class Template:
    """Class for storing Templates and associated information"""
    residues: list[Residue]
    pdb_id: str = field(default=None)
    mcsa_id: int = field(default=None)
    cluster: str = field(default=None)
    uniprot_id: str = field(default=None)
    organism: str = field(default=None)
    organism_id: Optional[int] = field(default=None)
    resolution: float = field(default=None)
    experimental_method: str = field(default=None)
    ec: str = field(default=None)

    @classmethod
    # TODO am I using self and cls correctly??
    # TODO should I use the parse functions like this - like static methods inside the class?
    def load(cls, file: TextIO, warn: bool = True) -> Template:
        atom_lines = []
        residues = []
        metadata = {}

        _PARSERS = {
            'PDB_ID': cls.parse_pdb_id,
            'UNIPROT_ID': cls.parse_uniprot_id,
            'MCSA_ID': cls.parse_mcsa_id,
            'CLUSTER': cls.parse_cluster,
            'ORGANISM_NAME': cls.parse_organism_name,
            'ORGANISM_ID': cls.parse_organism_id,
            'RESOLUTION': cls.parse_resolution,
            'EXPERIMENTAL_METHOD': cls.parse_experimental_method,
            'EC': cls.parse_ec
        }
                
        for line in filter(str.strip, file):
            tokens = line.split()
            if tokens[0] == 'REMARK':
                parser = _PARSERS.get(tokens[1])
                if parser is not None:
                    parser(tokens, metadata, warn=warn) # TODO parser functions do not need call to self since they fill metadata?
            elif tokens[0] == 'ATOM':
                atom_lines.append(line)
            else:
                continue

        for atom_triplet in chunks(atom_lines, 3): # yield chunks of 3 atom lines each
            atoms = tuple(Atom.loads(line) for line in atom_triplet)
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = { (atom.aa_type, atom.chain, atom.pdb_residue) for atom in atoms }
            if len(unique_residues) != 1:
                raise ValueError('Multiple residues found in atom lines')
            residues.append(Residue(atoms))

        return Template(
            residues = residues,
            **metadata, # unpack everything we parsed into metadata
        )

    @property
    def true_size(self) -> int: # Number of residues in the template as counted by triplets of atoms
        #TODO add int staements like this to other properties for auto-generating documentation
        """`int`: The number of residues in the template, as counted by triplets of atoms.
        """
        return len(self.residues)
    
    @cached_property
    def size(self) -> int: # Number of residues as evaluated, the effective size
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
        for residue in residues:
            if residue.specificity < 100 and residue.backbone == False: # type specific and not Backbone
                effective_size += 1
        return effective_size

    @cached_property
    def multimeric(self) -> bool: # if the template is split across multiple protein chains
        return all(res.chain == residues[0].chain for res in residues)

    @cached_property
    def relative_order(self) -> list: # list with length of deduplicated true_size
        # Now extract relative template order
        pdb_residue_list = [res.pdb_residue for res in residues]
        deduplicated_residue_list = list(dict.fromkeys(pdb_residue_list)) # remove duplicates while preserving order
        sorted_residues = sorted(deduplicated_residue_list)
        return [sorted_residues.index(i) for i in deduplicated_residue_list] # determine the relative order of residues

    def parse_pdb_id(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        metadata['pdb_id'] = tokens[2]

    def parse_uniprot_id(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', tokens[2])
        if match:
            metadata['uniprot_id'] = match.group()
        else:
            if warn:
                warnings.warn(f'Did not find a valid UniProt ID, found {tokens[2]}')

    def parse_mcsa_id(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['mcsa_id'] = int(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Did not find a M-CSA ID, found {tokens[2]}')

    def parse_cluster(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        metadata['cluster'] = tokens[2]

    def parse_organism_name(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        metadata['organism'] = tokens[2]

    def parse_organism_id(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['organism_id'] = int(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Ill-formatted organism ID: {tokens[2]}')

    def parse_resolution(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        try:
            metadata['resolution'] = float(tokens[2])
        except ValueError:
            if warn:
                warnings.warn(f'Ill-formatted pdb resolution: {tokens[2]}')

    def parse_experimental_method(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        metadata['experimental_method'] = tokens[2]

    def parse_ec(tokens: list[str], metadata: dict[str, object], warn: bool = True):
        match = re.search(r'\d{1,2}(\.(\-|\d{1,2})){3}', tokens[2])
        if match:
            metadata['ec'] = match.group()
        else:
            if warn:
                warnings.warn(f'Did not find a valid EC number, found {tokens[2]}')
                

def load_templates(template_dir=files(__package__).joinpath('jess_templates_20230210')) -> Iterator[Template]:
    """Load templates from a given directory, recursively.
    """
    if not Path(template_dir).exists():
        raise FileNotFoundError(template_dir)
    elif not Path(template_dir).is_dir():
        raise NotADirectoryError(template_dir)

    template_files = get_template_files(template_dir)
    for template_file in template_files:
        with template_file.open() as f:
            try:
                # maybe we want to yield the filepath to the template too? like as a tuple or something??
                yield Template.load(file=f)
            except ValueError as exc:
                raise ValueError(f"Failed to parse {template_file!r}") # from exc

if __name__ == "__main__":
    # TODO this is here for debugging purposes - delete later
    templates = list(load_templates())
    print(len(templates))
