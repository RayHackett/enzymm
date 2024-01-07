"""
written by Ray
script: dssp_run.py
this is called by batch_dssp.py: 
script to run dssp against a list of input structues and print the relative solvent accessibility to the occupancy column in the output pdb structure

needs Anaconda environment with Biopython and DSSP

#conda install -c salilab dssp
#conda install -c conda-forge "biopython>=1.80"  modified by Yannis!
#conda install -c anaconda pandas
#conda install -c anaconda requests 

Comments on Biopython - Version modified by Yannis:
# changes to DSSP.py see https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
# if chainid (line[11] == '>' then look at the end of the line)

Inputs:
    list of pdb structure files supplied with -i
    
Outputs:
    ./dssp_batch?.list with all the successfully dssp scored structures

"""

from Bio.PDB import Select, PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP
# see https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
from Bio.PDB.PDBIO import PDBIO
from pathlib import Path
import argparse
import traceback
import json
import warnings

# somehow this doesnt quite work ... dont know why but doesnt really matter
warnings.filterwarnings('ignore', '.*QCPSuperimposer module will be removed.*', )

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='file with list of protein structures to scan')
args = parser.parse_args()

pdb_file = args.input

cwd = Path.cwd()

results_path = Path(cwd,'dssp_results') 
results_path.mkdir(exist_ok=True)

pdb_path = Path(cwd,'cifchains') 
pdb_path.mkdir(exist_ok=True)

with open('/homes/hackett/catalytic_residues_homologues.json', 'r') as f:
    homologs = json.load(f)

pdb_chain_dict = {}
for i in homologs:
    for res in i['residue_chains']:
        #if res['assembly'] != None: # we ignore the ones with assembly == none
        name = res['pdb_id']
        if name not in pdb_chain_dict:
            pdb_chain_dict[name] = set()
        if res['is_reference']:
            pdb_chain_dict[name].add(res['assembly_chain_name'])
        else:
            pdb_chain_dict[name].add(res['chain_name'])

with open(pdb_file, 'r') as f:
    lines = f.readlines()
pdb_list = []
for line in lines:
    pdb_list.append(line.strip())
    
class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0

dssp_successes = []
for pdb in pdb_list:
    dssp_path = '/hps/software/users/thornton/hackett/anaconda3/envs/mod_biopy/bin/mkdssp'
    try:
        #structure = PDBParser(QUIET=True).get_structure(str(Path(pdb).name), Path(pdb))
        structure = MMCIFParser(QUIET=True).get_structure(str(Path(pdb).name), Path(pdb))
        # we get the correct chain from the homologs dictionary
        chains = list(pdb_chain_dict[str(Path(pdb).stem)[:4].lower()])
    
        model = structure[0] # always only one model per pdb file since we are dealing with predicted structures
        dssp = DSSP(model, Path(pdb), dssp=dssp_path, acc_array='Miller') # Miller since the Group useses that
        dssp_dict = dict(dssp)
        # keys are (chainid, ('icode', resseq, ' '))
        # However icode is returned empty for some reason. so (' ', 'resseq', ' ')
        # value is a tuple of:
        # (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
        # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
        # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    
        # add the relative solvent accessibility (normalized with Miller) to the occupancy column in the structure
        for i, residue in enumerate(structure.get_residues()):
                # check if it is a value
                data = dssp_dict.get((residue.get_full_id()[2], (' ', residue.id[1], residue.id[2])))
                if data and isinstance(data[3], float): # incomplete residues are not evaluated, HETATM get 'NA'
                    rsa = float(data[3])
                    #ss = data[2] # secondary structure
                else:
                    rsa = 9.99
                #plddt = residue['CA'].get_bfactor() / 100.0 # plddt
                atoms = residue.get_atoms()
                for atom in atoms:
                    atom.set_occupancy(rsa)
                    
        # save the modified structure to the same file location
        for chain in chains:
            pdb_chain_file = str(Path(pdb_path,str(Path(pdb).stem[:4].upper() + chain + '.pdb')))
            io=PDBIO(for_jess=True)
            io.set_structure(structure)
            io.save(pdb_chain_file, ChainSelect(chain), preserve_atom_numbering=True)
            with open(pdb_chain_file, 'r') as f:
                content = f.read()
            if content.strip() == 'END':
                Path(pdb_chain_file).unlink()
            else:
                dssp_successes.append(pdb_chain_file)
        
    except Exception:
        traceback.print_exc()
            
with open(Path(results_path, Path(pdb_file).name.split('.')[0] + '.list'), 'w') as f:
    for i in dssp_successes:
        f.write(str(i) + '\n')
    
