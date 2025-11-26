from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List, Dict
import warnings
import sys
import os

import rich.console
from pyjess import Molecule

from enzymm import __version__
from enzymm.template import load_templates
from enzymm.jess_run import Matcher, Match, load_molecules


def build_parser() -> argparse.ArgumentParser:
    """Parse Arguments with Argparse. Returns args object"""

    class ReadListAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            with values.open("r") as f:  # type: ignore
                for line in f:
                    dest = getattr(namespace, self.dest)
                    dest.append(Path(line.strip()))

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage="""
            Minimal use: enzymm -i query.pdb -o result.tsv
            
            Recommended use: enzymm -l query_pdbs.list -o results.tsv -v --pdbs pdb_folder --include-query

            """,
        description=(
            f"EnzyMM - The Enzyme Motif Miner - version {__version__}\n\n"
            "Geometric matching of catalytic motifs in protein structures \n\n"
            "MIT License\n\n"
            "Copyright (c) 2025 Raymund Hackett <r.e.hackett@lumc.nl>\n"
        ),
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"enzymm {__version__}"
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=(
            len(os.sched_getaffinity(0))
            if sys.platform == "linux"
            else os.cpu_count() or 1
        ),  # len(os.sched_getaffinity(0)) doesnt work on mac/win
        help="The number of threads to spawn for jobs in parallel. Pass 0 to select all cores. Negative numbers: leave this many cores free.",
    )

    group = parser.add_argument_group("Mandatory Parameters")
    group.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output tsv file to which results should get written",
    )

    group = parser.add_argument_group("Inputs - Use either or combine")
    # inputs: either a list of paths, or directly a path (or any combination)
    group.add_argument(
        "-i",
        "--input",
        type=Path,
        help="File path to a PDB or mmCIF file to use as query",
        action="append",
        dest="files",
        default=[],
    )
    group.add_argument(
        "-l",
        "--list",
        type=Path,
        help="File containing a list of PDB or mmCIF files to read",
        action=ReadListAction,
        dest="files",
    )

    group = parser.add_argument_group("Optional Arguments")
    # optional arguments
    group.add_argument(
        "--pdbs",
        type=Path,
        help="Output directory to which results should get written",
        default=None,
    )
    group.add_argument(
        "-p",
        "--parameters",
        nargs=3,
        default=None,
        type=float,
        help="Fixed Jess parameters for all templates. Jess space seperated parameters rmsd, distance, max_dynamic_distance",
    )
    group.add_argument(
        "-t",
        "--template-dir",
        type=Path,
        default=None,
        help="Path to directory containing jess templates. This directory will be recursively searched.",
    )
    group.add_argument(
        "-c",
        "--conservation-cutoff",
        type=float,
        default=0,
        help="Atoms with a value in the B-factor column below this cutoff will be excluded form matching to the templates. Useful for predicted structures.",
    )

    group = parser.add_argument_group("Flags")
    group.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="If process information and time progress should be printed to the command line",
    )
    group.add_argument(
        "-w",
        "--warn",
        default=False,
        action="store_true",
        help="If warings about bad template processing or suspicous and missing annotations should be raised",
    )
    group.add_argument(
        "-q",
        "--include-query",
        default=False,
        action="store_true",
        help="Include the query structure together with the hits in the pdb output",
    )
    group.add_argument(
        "--include-template",
        default=False,
        action="store_true",
        help="Include the template structure together with the hits in the pdb output",
    )
    group.add_argument(
        "-u",
        "--unfiltered",
        default=False,
        action="store_true",
        help="If set, matches which logistic regression predicts as false based on RMSD and resdiue orientation will be retained. By default, matches predicted as false are removed.",
    )
    group.add_argument(
        "--transform",
        default=False,
        action="store_true",
        help="If set, one pdb file per matched template pdb with will be written in the coordinate system of that template",
    )
    group.add_argument(
        "--skip-smaller-hits",
        default=False,
        action="store_true",
        help="If set, will not search with smaller templates if larger templates have already found hits.",
    )
    group.add_argument(
        "--match-small-templates",
        default=False,
        action="store_true",
        help="If set, templates with less then 3 defined sidechain residues will still be matched.",
    )
    group.add_argument(
        "--skip-annotation",
        default=False,
        action="store_true",
        help="If set, M-CSA derived templates will NOT be annotated with extra information.",
    )
    return parser


def main(argv: Optional[List[str]] = None, stderr=sys.stderr):

    parser = build_parser()
    args = parser.parse_args(args=argv)
    if not args.files:
        raise ValueError("No input files were passed. Use -i and/or -l.")

    jess_params = None
    if args.parameters:
        jess_params_list = [i for i in args.parameters]
        # jess parameters
        # we use different parameters for different template residue numbers - higher number more generous parameters
        rmsd = jess_params_list[0]  # in Angstrom, typcically set to 2
        distance = jess_params_list[
            1
        ]  # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
        max_dynamic_distance = jess_params_list[
            2
        ]  # if equal to distance dynamic is off: this option is currenlty dysfunctional

        known_distances = [float(i) for i in Match.ensemble_model.ensemble[3].keys()]
        if distance not in known_distances and not args.unfiltered:
            raise ValueError(
                "Filtering paraterms only established for pairwise distances in 0.7, 0.8, ... , 2.0A "
            )

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
        molecules = load_molecules(
            molecule_paths=args.files,
            conservation_cutoff=args.conservation_cutoff,
        )

        templates = list(
            load_templates(
                template_dir=args.template_dir,
                warn=args.warn,
                verbose=args.verbose,
                with_annotations=not args.skip_annotation,
            )
        )

        ############ Initialize Matcher object ################################
        matcher = Matcher(
            templates=templates,
            jess_params=jess_params,
            warn=args.warn,
            verbose=args.verbose,
            skip_smaller_hits=args.skip_smaller_hits,
            match_small_templates=args.match_small_templates,
            cpus=args.jobs,
            filter_matches=not args.unfiltered,
            console=rich.console.Console(file=stderr),
        )

        ############ Call Matcher.run ##########################################
        processed_molecules = matcher.run(molecules=molecules)

        ######### Writing Output ##########################################
        out_tsv = args.output

        if not out_tsv.parent.exists():
            out_tsv.parent.mkdir(parents=True, exist_ok=True)
            if args.warn:
                warnings.warn(f"{out_tsv.parent.resolve()} dir to output was created")
        elif out_tsv.exists() and args.warn:
            warnings.warn(f"The output file {out_tsv.resolve()} will be overwritten!")

        if args.verbose:
            print(f"Writing output to {out_tsv.resolve()}")
            print(
                f"Matches predicted by logistic regression as false are {'' if args.unfiltered else 'not '}reported"
            )

        with open(out_tsv, "w", newline="", encoding="utf-8") as tsvfile:
            tsvfile.write(f"# Command: {' '.join(sys.argv)}\n")
            for index, (molecule, matches) in enumerate(processed_molecules.items()):
                for jndex, match in enumerate(matches):
                    i = index + jndex
                    match.index = jndex + 1  # 1 indexed matches per query
                    match.dump(
                        tsvfile,
                        header=(i == 0),
                    )  # one line per match, write header only for the first match

        def write_hits2pdb(matches: List[Match], filename: str, outdir: Path):
            # make sure molecule().id is unique!
            # TODO fix model indx upon multiple calls
            with open(
                Path(outdir, f"{filename}_matches.pdb"), "a", encoding="utf-8"
            ) as pdbfile:
                if args.include_query:
                    # write the molecule structure to the top of the pdb output too
                    pdbfile.write("MODEL        0\n")
                    matches[0].dump_query(file=pdbfile, transform=args.transform)
                    pdbfile.write("ENDMDL\n\n")

                model_idx = 1
                for match in matches:
                    pdbfile.write(f"MODEL        {model_idx}\n")
                    match.dump2pdb(pdbfile, transform=args.transform)
                    model_idx += 1
                    pdbfile.write("\n")

                if args.include_template:
                    for match in matches:
                        pdbfile.write(f"MODEL        {model_idx}\n")
                        match.dump_template(pdbfile, transform=args.transform)
                        pdbfile.write("ENDMDL\n\n")
                        model_idx += 1

        if args.pdbs:
            args.pdbs.mkdir(parents=True, exist_ok=False)

            if args.transform:
                t_dict: Dict[str, Dict[Molecule, List[Match]]] = {}
                for molecule, matches in processed_molecules.items():
                    for match in matches:
                        pdb_id = match.hit.template().pdb_id
                        if pdb_id not in t_dict:
                            t_dict[pdb_id] = {}
                        if molecule not in t_dict[pdb_id]:
                            t_dict[pdb_id][molecule] = []
                        t_dict[pdb_id][molecule].append(match)

                for template_pdb, query_dict in t_dict.items():
                    for query_mol, matches in query_dict.items():
                        write_hits2pdb(
                            matches=matches, filename=template_pdb, outdir=args.pdbs
                        )

            else:
                # everything in the query reference frame
                # write one pdb per query
                for molecule, matches in processed_molecules.items():
                    write_hits2pdb(
                        matches=matches,
                        filename=molecule.id,  # type: ignore
                        outdir=args.pdbs,
                    )

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
