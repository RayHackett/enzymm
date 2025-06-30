from matcher.template import load_templates
from pathlib import Path
from typing import Set, List


def get_list_of_template_pdbchains() -> List[str]:

    templates = list(load_templates(template_dir=None, with_annotations=False))

    template_pdb_chains: Set[str] = set()
    # template_ec_numbers: Set[str] = set()
    # template_mcsa_numbers: Set[str] = set()

    for template in templates:
        # if template.effective_size >= 3:
        template_resids = template.residues
        chains: Set[str] = set()
        for resid in template_resids:
            chains.add(resid.chain_id)

        for chain in chains:
            template_pdb_chains.add(template.pdb_id.lower() + "." + chain)  # type: ignore

        # template_ec_numbers.update(template.ec)
        # template_mcsa_numbers.add(template.mcsa_id) # type: ignore

    return template_pdb_chains


def get_list_of_template_ecs() -> List[str]:
    templates = list(load_templates(template_dir=None, with_annotations=False))

    template_ec_numbers: Set[str] = set()

    for template in templates:
        template_ec_numbers.update(template.ec)

    return template_ec_numbers


def write_list_to_file(list: List[str], filename: Path):

    with open(filename, "w") as f:
        for i in list:
            f.write(i + "\n")


if __name__ == "__main__":
    template_pdbchains = get_list_of_template_pdbchains()
    write_list_to_file(
        list=template_pdbchains, filename=Path("template_pdbchains.list")
    )
