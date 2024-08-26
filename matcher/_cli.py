from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List, Dict, TYPE_CHECKING
from collections import defaultdict
import warnings
import sys
import pyjess
from .template import load_templates
from . import __version__

if TYPE_CHECKING:
    from .jess_run import Matcher, Match


def build_parser() -> argparse.ArgumentParser:
    """Parse Arguments with Argparse. Returns args object"""

    class ReadListAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            with values.open("r") as f:
                for line in f:
                    dest = getattr(namespace, self.dest)
                    dest.append(Path(line.strip()))

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, exit_on_error=False
    )
    # positional arguments
    # parser.add_argument('input', type=Path, help='File path of a single query pdb file OR File path of a file listing multiple query pdb files seperated by linebreaks')
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output tsv file to which results should get written",
    )
    parser.add_argument(
        "--pdbs",
        type=Path,
        help="Output directory to which results should get written",
        default=None,
    )

    # inputs: either a list of paths, or directly a path (or any combination)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        help="File path to a PDB file to use as query",
        action="append",
        dest="files",
        default=[],
    )
    parser.add_argument(
        "-l",
        "--list",
        type=Path,
        help="File containing a list of PDB files to read",
        action=ReadListAction,
        dest="files",
    )

    # optional arguments
    parser.add_argument(
        "-j",
        "--jess",
        nargs=3,
        default=None,
        type=float,
        help="Fixed Jess parameters for all templates. Jess space seperated parameters rmsd, distance, max_dynamic_distance",
    )
    parser.add_argument(
        "-t",
        "--template-dir",
        type=Path,
        default=None,
        help="Path to directory containing jess templates. This directory will be recursively searched.",
    )
    parser.add_argument(
        "-c",
        "--conservation-cutoff",
        type=float,
        default=0,
        help="Atoms with a value in the B-factor column below this cutoff will be excluded form matching to the templates. Useful for predicted structures.",
    )
    parser.add_argument(
        "-n",
        "--n-jobs",
        type=int,
        default=0,
        help="The number of threads to run in parallel. Pass 0 to select all. Negative numbers: leave this many threads free.",
    )

    # Boolean optional arguments
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="If process information and time progress should be printed to the command line",
    )
    parser.add_argument(
        "-w",
        "--warn",
        default=False,
        action="store_true",
        help="If warings about bad template processing or suspicous and missing annotations should be raised",
    )
    parser.add_argument(
        "-q",
        "--include-query",
        default=False,
        action="store_true",
        help="Include the query structure together with the hits in the pdb output",
    )
    parser.add_argument(
        "-f",
        "--filter-matches",
        default=True,
        action="store_false",
        help="If set, matches which logistic regression predicts as false based on RMSD and resdiue orientation will be retained. By default, matches predicted as false are removed",
    )
    parser.add_argument(
        "--transform",
        default=False,
        action="store_true",
        help="If set, will transform the coordinate system of the hits to that of the template in the pdb output",
    )
    parser.add_argument(
        "--skip-smaller-hits",
        default=False,
        action="store_true",
        help="If set, will not search with smaller templates if larger templates have already found hits.",
    )
    parser.add_argument(
        "--match-small-templates",
        default=False,
        action="store_true",
        help="If set, templates with less then 3 defined sidechain residues will still be matched.",
    )
    return parser


def main(argv: Optional[List[str]] = None, stderr=sys.stderr):

    parser = build_parser()
    args = parser.parse_args(args=argv)
    if not args.files:
        parser.error(
            "No input files were passed. Use the -i and/or -l flags to pass input."
        )

    jess_params = None
    if args.jess:
        jess_params_list = [i for i in args.jess]
        # jess parameters
        # we use different parameters for different template residue numbers - higher number more generous parameters
        rmsd = jess_params_list[0]  # in Angstrom, typcically set to 2
        distance = jess_params_list[
            1
        ]  # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
        max_dynamic_distance = jess_params_list[
            2
        ]  # if equal to distance dynamic is off: this option is currenlty dysfunctional

        jess_params = {
            3: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
            4: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
            5: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
            6: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
            7: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
            8: {
                "rmsd": rmsd,
                "distance": distance,
                "max_dynamic_distance": max_dynamic_distance,
            },
        }

    try:

        ######## Loading all query molecules #############################
        molecules: List[pyjess.Molecule] = []
        stem_counter: Dict[str, int] = defaultdict(int)
        for molecule_path in args.files:
            stem = Path(molecule_path).stem
            stem_counter[stem] += 1
            if stem_counter[stem] > 1:
                # In case the same stem occurs multiple times, create a unique ID using the stem and a running number starting from 2
                unique_id = f"{stem}_{stem_counter[stem]}"
            else:
                unique_id = stem

            molecule = pyjess.Molecule.load(
                str(molecule_path), id=unique_id
            )  # by default it will stop at ENDMDL
            if args.conservation_cutoff:
                molecule.conserved(args.conservation_cutoff)
                # molecule = molecule.conserved(args.conservation_cutoff)
                # conserved is a method called on a molecule object that returns a filtered molecule
                # atoms with a temperature-factor BELOW the conservation cutoff will be excluded
            if molecule:
                molecules.append(
                    molecule
                )  # load a molecule and filter it by conservation_cutoff
            elif args.warn:
                warnings.warn(f"received an empty molecule from {molecule_path}")

        if not molecules and args.warn:
            warnings.warn("received no molecules from input")

        if args.template_dir:
            templates = list(
                load_templates(
                    template_dir=Path(args.template_dir),
                    warn=args.warn,
                    verbose=args.verbose,
                    cpus=args.n_jobs,
                )
            )
        else:
            templates = list(
                load_templates(warn=args.warn, verbose=args.verbose, cpus=args.n_jobs)
            )

        ############ Initialize Matcher object ################################
        matcher = Matcher(
            templates=templates,
            jess_params=jess_params,
            conservation_cutoff=args.conservation_cutoff,
            warn=args.warn,
            verbose=args.verbose,
            skip_smaller_hits=args.skip_smaller_hits,
            match_small_templates=args.match_small_templates,
            cpus=args.n_jobs,
            filter_matches=args.filter_matches,
        )

        ############ Call Matcher.run ##########################################
        processed_molecules = matcher.run(molecules=molecules)

        ######### Writing Output ##########################################
        out_tsv = args.output

        if not out_tsv.parent.exists():
            out_tsv.parent.mkdir(parents=True, exist_ok=True)
            if args.warn:
                warnings.warn(
                    f"{out_tsv.parent.resolve()} directory tree to output tsv file did not exist and was created"
                )
        elif out_tsv.exists() and args.warn:
            warnings.warn(
                f"The specified output tsv file {out_tsv.resolve()} already exists and will be overwritten!"
            )

        if args.verbose:
            print(f"Writing output to {out_tsv.resolve()}")
            print(
                f"Matches predicted by logistic regression as false are {'not ' if args.filter_matches else ''}reported"
            )

        with open(out_tsv, "w", newline="", encoding="utf-8") as tsvfile:
            tsvfile.write(f"# Version {__version__}")
            for index, (molecule, matches) in enumerate(processed_molecules.items()):
                for jndex, match in enumerate(matches):
                    i = index + jndex
                    match.index = jndex + 1  # 1 indexed matches per query
                    match.dump(
                        tsvfile, header=i == 0
                    )  # one line per match, write header only for the first match

        def write_hits2pdb(matches: List[Match], filename: str, outdir: Path):
            outdir.mkdir(parents=True, exist_ok=True)
            # make sure molecule.id is unique!
            with open(
                Path(outdir, f"{filename}_matches.pdb"), "w", encoding="utf-8"
            ) as pdbfile:
                if (
                    args.include_query
                ):  # write the molecule structure to the top of the pdb output too
                    for i, match in enumerate(matches):
                        match.dump2pdb(
                            pdbfile, include_query=i == 0, transform=args.transform
                        )
                else:
                    for match in matches:
                        match.dump2pdb(pdbfile, transform=args.transform)

        if args.pdbs:
            for molecule, matches in processed_molecules.items():
                write_hits2pdb(matches=matches, filename=molecule.id, outdir=args.pdbs)  # type: ignore

    except IsADirectoryError as exc:
        print("File is a directory:", exc.filename, file=stderr)
        return exc.errno

    except FileNotFoundError as exc:
        print("Failed to find file:", exc.filename, file=stderr)
        return exc.errno

    except FileExistsError as exc:
        print("File already exists:", exc.filename, file=stderr)
        return exc.errno

    return 0


if __name__ == "__main__":
    sys.exit(main())
