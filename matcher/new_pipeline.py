from pathlib import Path
from typing import List, Set, Iterator
from tqdm import tqdm

from .jess_run import matcher_run

def run_pipeline(list_of_queries_file: Path,  rmsd: float, distance: float, max_dynamic_distance: float, score_cutoff: float, outdir: Path, complete_search: bool):

    def read_line_by_line(file: Path) -> Iterator[Path]:
        try:
            for path in open(list_of_queries_file):
                yield path
        except FileNotFoundError:
            print(f"Input file not found: {list_of_queries_file}")

    if not outdir.is_dir():
        # TODO would it be better to create the directory instead and print a statement that a directory was created?
        raise FileNotFoundError(f"'{outdir}' is not an existing directory.")
        
    # Pass all the molecule paths as a list to the matcher_run function
    for molecule_path in tqdm(read_line_by_line(list_of_queries_file), ncols=100, description='Searching Query structures against Template library'):
        if molecule_path.exists():
            # check if its a .pdb or .ent file
            if molecule_path.suffix.lower() not in {".pdb", ".ent"}:
                raise ValueError(f'File {path} did not have expected .pdb or .ent extensions')

            matcher_run(query_file=query_file,  rmsd=rmsd, distance=distance, max_dynamic_distance=max_dynamic_distance, score_cutoff=score_cutoff, outdir=outdir, complete_search=complete_search)
            # TODO should I hvae matcher_run return something? true or fasle based on hit found or not?
            # should matcher_run return a Path to the file it wrote?

            # TODO how to handle queries for which no hits were found?
            # with open(new_working_file, 'w') as f:
            #     for i in queries_with_no_hits:
            #         f.write(i + '\n')

            # TODO aggregate all the tsv files that were written into one big one?

        else:
            raise FileNotFoundError(f'Molecule file {molecule_path} does not exist')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, help='file containing a list of pdb files to search seperated by linebreaks.')
    parser.add_argument('-j', '--jess', nargs = '+', help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff, optional flags as a string')
    parser.add_argument('-o', '--output', type=Path, help='Output Directory to which results should get written', default=Path.cwd())
    parser.add_argument('-c', '--complete', tuype=bool, help='If True continue search with smaller templates even after hits with larger template have been found', default=True)
    args = parser.parse_args()
    
    list_of_queries_file = args.input
    jess_params = [i for i in args.jess]
    outdir = args.output
    complete_search = args.complete

    # jess parameters
    # we use different parameters for different template residue numbers - higher number more generous parameters
    rmsd = jess_params[0] # in Angstrom, typcically set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    score_cutoff = jess_params[3] # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example
    
    if len(jess_params) != 4:
        sys.exit('Wrong number of Jess Parameters') # TODO what kind of error should I throw?

    run_pipeline(list_of_queries_file, rmsd, distance, max_dynamic_distance, score_cutoff, outdir, complete_search)