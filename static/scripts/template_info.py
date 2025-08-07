import sys
import argparse
from pathlib import Path
from typing import Set, List, Optional

from enzymm.template import load_templates


def build_parser() -> argparse.ArgumentParser:
    """Parse Arguments with Argparse. Returns args object"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=""" python template_info -t pdbids -d comma
            """,
        description="Write information about templates to a file",
    )
    parser.add_argument(
        "-t",
        "--type",
        help="Which type of information should be collected",
        choices=["pdbids", "pdbid-assembly", "pdbchains", "ec"],
        required=True,
    )

    parser.add_argument(
        "-d",
        "--deliminator",
        help="Which deliminator should be used",
        choices=["newline", "tab", "comma", "space"],
        required=False,
        default="comma",
    )

    return parser


def get_list_of_template_pdbids() -> List[str]:

    templates = list(load_templates(template_dir=None, with_annotations=False))

    template_pdbids: Set[str] = set()

    for template in templates:
        template_pdbids.add(template.pdb_id.lower())  # type: ignore

    return list(template_pdbids)


def get_list_of_template_pdbids_with_assembly() -> List[str]:

    templates = list(load_templates(template_dir=None, with_annotations=True))

    template_pdbids: Set[str] = set()

    for template in templates:

        template_pdbids.add(template.pdb_id.lower() + "-" + str(template.assembly))  # type: ignore

    return list(template_pdbids)


def get_list_of_template_pdbchains() -> List[str]:

    templates = list(load_templates(template_dir=None, with_annotations=False))

    template_pdb_chains: Set[str] = set()

    for template in templates:
        # if template.effective_size >= 3:
        template_resids = template.residues
        chains: Set[str] = set()
        for resid in template_resids:
            chains.add(resid.chain_id)

        for chain in chains:
            template_pdb_chains.add(template.pdb_id.lower() + "." + chain)  # type: ignore

    return list(template_pdb_chains)


def get_list_of_template_ecs() -> List[str]:
    templates = list(load_templates(template_dir=None, with_annotations=False))

    template_ec_numbers: Set[str] = set()

    for template in templates:
        template_ec_numbers.update(template.ec)

    return list(template_ec_numbers)


def write_list_to_file(
    l: List[str],
    filename: Path,
    deliminator: str,
):

    if deliminator == "newline":
        del_char = "\n"
    elif deliminator == "tab":
        del_char = "\t"
    elif deliminator == "comma":
        del_char = ","
    elif deliminator == "space":
        del_char = " "

    with open(filename, "w") as f:
        for i in l:
            f.write(i)
            f.write(del_char)


def main(argv: Optional[List[str]] = None, stderr=sys.stderr):

    parser = build_parser()
    args = parser.parse_args(args=argv)

    if args.type == "pdbids":
        item_list = get_list_of_template_pdbids()
    elif args.type == "pdbchains":
        item_list = get_list_of_template_pdbchains()
    elif args.type == "pdbid-assembly":
        item_list = get_list_of_template_pdbids_with_assembly()
    elif args.type == "ec":
        item_list = get_list_of_template_ecs()
    else:
        raise ValueError(f"Unknown type {args.type}")

    name = args.type
    write_list_to_file(
        l=item_list,
        filename=Path(f"template_{name}.list"),
        deliminator=args.deliminator,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
