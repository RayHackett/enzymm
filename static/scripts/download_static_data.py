from pathlib import Path
from typing import List
import json
import pandas as pd  # type: ignore
import urllib.request
import gzip
import shutil
from collections import OrderedDict
from dataclasses import asdict
import requests  # type: ignore
from requests.adapters import HTTPAdapter  # type: ignore
from urllib3.util.retry import Retry  # type: ignore
import time

from enzymm.utils import SetEncoder  # type: ignore
from template_info import get_list_of_template_pdbchains

from enzymm.mcsa_info import (
    HomologousPDB,
    ReferenceCatalyticResidue,
    NonReferenceCatalyticResidue,
)


def request_url(
    url: str, acceptable_stati: List[int], timeout: int = 10, max_retries: int = 10
):
    """
    Fetch content from a url, while handling retries and url issues.

    Arguments:
        acceptable_stati: List of HTTP Status codes to accept. If Code is is accepted but not 200, return an empty string.
        timeout: time to sleep between retries
        max_retries: number of times to retry reaching a url

    Returns:
        `Response` object or empty string or prints a warning.
    """
    # TODO figure out maxretries

    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    r = session.get(url)
    request_counter = 0
    while request_counter < max_retries and r.status_code in range(
        500, 600
    ):  # server errors 500-599; 429 too many requests
        time.sleep(timeout)
        r = session.get(url)
        request_counter += 1

    if r.status_code == 429:
        time.sleep(int(r.headers["Retry-After"]))

    if r.status_code in acceptable_stati:
        if r.status_code == 200:
            return r
        else:
            return ""
    else:
        print("Failed to get url", url, " with status code: ", r.status_code)


def check_other_files(outdir: Path):
    Path(outdir).mkdir(parents=True, exist_ok=True)

    EC_cofactor_mapping = Path(outdir, "EClist_cofactors_forRH.csv")
    # AlphaFold_CATH = Path(outdir, 'cath-v4_3_0.alphafold-v2.2022-11-22.tsv')
    PDBchain_to_CATH_Uniprot = Path(outdir, "pdb_chain_cath_uniprot.csv")
    PDBchain_to_EC_Uniprot = Path(outdir, "pdb_chain_enzyme.csv")
    # PDBchain_to_InterPro = Path(outdir, "pdb_chain_interpro.csv")
    CATH_names_mapping = Path(outdir, "cath-domain-list.txt")

    if not Path(EC_cofactor_mapping).is_file():
        raise FileNotFoundError(
            f"EC to cofactor mapping csv expected at {EC_cofactor_mapping.resolve()} is missing!"
        )

    if not Path(PDBchain_to_CATH_Uniprot).is_file():
        compressed = Path(outdir, "pdb_chain_cath_uniprot.csv.gz")
        urllib.request.urlretrieve(
            "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_cath_uniprot.csv.gz",
            compressed,
        )

        with gzip.open(compressed, "rb") as f_gz:
            with open(PDBchain_to_CATH_Uniprot, "wb") as f:
                shutil.copyfileobj(f_gz, f)
        # Mapping PDBchains to CATH name and Uniprot! Obtain via:
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

        Path.unlink(compressed)

    if not Path(PDBchain_to_EC_Uniprot).is_file():
        compressed = Path(outdir, "pdb_chain_enzyme.csv.gz")

        urllib.request.urlretrieve(
            "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_enzyme.csv.gz",
            compressed,
        )

        with gzip.open(compressed, "rb") as f_gz:
            with open(PDBchain_to_EC_Uniprot, "wb") as f:
                shutil.copyfileobj(f_gz, f)

        Path.unlink(compressed)

        # Mapping PDBchains to EC and Uniprot
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

    # if not Path(PDBchain_to_InterPro).is_file():
    #     compressed = Path(outdir, "pdb_chain_interpro.csv.gz")

    #     urllib.request.urlretrieve(
    #         "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_interpro.csv.gz",
    #         compressed,
    #     )

    #     with gzip.open(compressed, "rb") as f_gz:
    #         with open(PDBchain_to_InterPro, "wb") as f:
    #             shutil.copyfileobj(f_gz, f)

    #     Path.unlink(compressed)

    #     # Mapping PDBchains to EC and Uniprot
    #     # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

    if not Path(CATH_names_mapping).is_file():
        # Mapping CATH ids to CATH numbers and CATH domain names
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt",
            CATH_names_mapping,
        )

    # if not Path(AlphaFold_CATH).is_file():
    #     # Mapping UniProt AFDB Proteins to CATH Domains
    #     urllib.request.urlretrieve(
    #       'ftp://orengoftp.biochem.ucl.ac.uk/alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv',
    #       AlphaFold_CATH
    #     )


def make_pdb_sifts_df(outdir: Path):

    dir = Path("/home/ray/Documents/template_matching/enzymm/data")

    PDBchain_to_CATH_Uniprot = Path(dir, "pdb_chain_cath_uniprot.csv")
    PDBchain_to_EC_Uniprot = Path(dir, "pdb_chain_enzyme.csv")
    # PDBchain_to_InterPro = Path(dir, "pdb_chain_interpro.csv")
    CATH_names_mapping = Path(dir, "cath-domain-list.txt")

    # load them all as df
    pdb_cath_uniprot_df = pd.read_csv(
        PDBchain_to_CATH_Uniprot, comment="#", dtype=str, na_filter=False
    )
    pdb_enzyme_uniprot_df = pd.read_csv(
        PDBchain_to_EC_Uniprot, comment="#", dtype=str, na_filter=False
    )
    # pdb_interpro_uniprot_df = pd.read_csv(
    #     PDBchain_to_InterPro, comment="#", dtype=str, na_filter=False
    # )
    cath_naming_df = pd.read_csv(
        CATH_names_mapping,
        comment="#",
        sep=r"\s+",
        header=None,
        dtype=str,
        na_filter=False,
    )

    # join the right columns to form the CATH NUMBER
    cath_naming_df["CATH_NUMBER"] = (
        cath_naming_df[1]
        + "."
        + cath_naming_df[2]
        + "."
        + cath_naming_df[3]
        + "."
        + cath_naming_df[4]
    )
    # rename a column
    cath_naming_df.rename(columns={0: "CATH_ID"}, inplace=True)
    pdb_cath_uniprot_df = pdb_cath_uniprot_df.merge(
        cath_naming_df[["CATH_ID", "CATH_NUMBER"]], on="CATH_ID"
    )  # add in the CATH number by CATH_ID

    pdb_enzyme_uniprot_df["PDBCHAIN"] = (
        pdb_enzyme_uniprot_df["PDB"] + pdb_enzyme_uniprot_df["CHAIN"]
    )
    pdb_cath_uniprot_df["PDBCHAIN"] = (
        pdb_cath_uniprot_df["PDB"] + pdb_cath_uniprot_df["CHAIN"]
    )
    # pdb_interpro_uniprot_df["PDBCHAIN"] = (
    #     pdb_interpro_uniprot_df["PDB"] + pdb_interpro_uniprot_df["CHAIN"]
    # )

    template_pdbchains = set(
        "".join(i.strip().split(".")) for i in get_list_of_template_pdbchains()
    )

    pdb_enzyme_uniprot_df = pdb_enzyme_uniprot_df[
        pdb_enzyme_uniprot_df["PDBCHAIN"].isin(template_pdbchains)
    ]
    pdb_cath_uniprot_df = pdb_cath_uniprot_df[
        pdb_cath_uniprot_df["PDBCHAIN"].isin(template_pdbchains)
    ]
    # pdb_interpro_uniprot_df = pdb_interpro_uniprot_df[pdb_interpro_uniprot_df['PDBCHAIN'].isin(template_pdbchains)]

    sifts_dict = {}
    for row in pdb_enzyme_uniprot_df.itertuples():
        if row.PDBCHAIN not in sifts_dict:
            sifts_dict[row.PDBCHAIN] = {
                "uniprot_id": row.ACCESSION,
                "cath": set(),
                "ec": set(),
                # "interpro": set(),
            }
        if row.EC_NUMBER != "?":
            sifts_dict[row.PDBCHAIN]["ec"].add(row.EC_NUMBER)
    for row in pdb_cath_uniprot_df.itertuples():
        if row.PDBCHAIN not in sifts_dict:
            sifts_dict[row.PDBCHAIN] = {
                "uniprot_id": row.SP_PRIMARY,
                "cath": set(),
                "ec": set(),
                # "interpro": set(),
            }
        sifts_dict[row.PDBCHAIN]["cath"].add(row.CATH_NUMBER)
    # for row in pdb_interpro_uniprot_df.itertuples():
    #     if row.PDBCHAIN not in sifts_dict:
    #         sifts_dict[row.PDBCHAIN] = {
    #             "uniprot_id": None,
    #             "cath": set(),
    #             "ec": set(),
    #             "interpro": set(),
    #         }
    #     sifts_dict[row.PDBCHAIN]["interpro"].add(row.INTERPRO_ID)

    # merge the dataframes
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with open(Path(outdir, "pdb_sifts.json"), "w") as f:
        json.dump(sifts_dict, f, indent=4, sort_keys=True, cls=SetEncoder)


def download_mcsa_entry_data(outdir: Path):
    file = Path(outdir, "mcsa_entry_info.json")
    if not file.is_file():
        url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page=1"
        next_exists = True
        entry_dict = OrderedDict()
        while next_exists:
            result = request_url(url, [200])
            general_page = result.json()
            for i in general_page["results"]:
                entry_dict[i["mcsa_id"]] = i
            if result.next:
                url = result.next
            else:
                next_exists = False

        Path(outdir).mkdir(parents=True, exist_ok=True)
        with open(file, "w") as f:
            json.dump(entry_dict, f, indent=4, sort_keys=True)


def get_mcsa_ec(outdir: Path):
    # m-csa entry map to EC number
    # every m-csa entry only has one EC number
    file = Path(outdir, "MCSA_EC_mapping.json")
    if not file.is_file():
        ec_dict = {}

        entry_file = Path(outdir, "mcsa_entry_info.json")
        if entry_file.is_file():
            with open(Path(outdir, "mcsa_entry_info.json"), "r") as f:
                entry_info = json.load(f)
        else:
            download_mcsa_entry_data(outdir=outdir)
            with open(Path(outdir, "mcsa_entry_info.json"), "r") as f:
                entry_info = json.load(f)

        for mcsa_id, results in entry_info.items():
            if mcsa_id not in ec_dict:
                ec_dict[mcsa_id] = ""
            if len(results["all_ecs"]) > 1:
                print(f"More than one EC found for Entry {mcsa_id}")
            for ec_num in results["all_ecs"]:
                ec_dict[mcsa_id] = ec_num

            Path(outdir).mkdir(parents=True, exist_ok=True)
            with open(file, "w") as f:
                json.dump(ec_dict, f, indent=4, sort_keys=True)


def get_mcsa_cath(outdir: Path):
    # m-csa entry map to CATH id
    # entry have multiple CATH ids
    file = Path(outdir, "MCSA_CATH_mapping.json")
    if not file.is_file():
        mcsa_cath = OrderedDict()

        entry_file = Path(outdir, "mcsa_entry_info.json")
        if entry_file.is_file():
            with open(Path(outdir, "mcsa_entry_info.json"), "r") as f:
                entry_info = json.load(f)
        else:
            download_mcsa_entry_data(outdir=outdir)
            with open(Path(outdir, "mcsa_entry_info.json"), "r") as f:
                entry_info = json.load(f)

        for mcsa_id, results in entry_info.items():
            cath_ids = set()
            for resid in results["residues"]:
                for chain in resid["residue_chains"]:
                    cath_ids.add(chain["domain_cath_id"])

            if "" in cath_ids:
                cath_ids.remove("")
            mcsa_cath[mcsa_id] = list(cath_ids)

        Path(outdir).mkdir(parents=True, exist_ok=True)
        with open(file, "w") as f:
            json.dump(mcsa_cath, f, indent=4, sort_keys=True)


def download_mcsa_homolog_info(outdir: Path):
    # download the information on catalytic residues found by homology to M-CSA entries
    file = Path(outdir, "catalytic_residue_homologs_information.json")

    # load the templates without annotations to get a list of pdbchains
    # NOTE that importing get_list_of_template_pdbchains consitutes a ciruclar dependency
    # ideally you therefore save this to a new file
    # while still keeping the original in the data folder
    template_pdbchains = set(
        "".join(i.strip().split(".")) for i in get_list_of_template_pdbchains()
    )

    # download the file from m-csa
    if not file.is_file():
        url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json"
        result = request_url(url, [200])
        homology_info = result.json()

        problematic_entries = set()
        ordered_homology_info: OrderedDict = OrderedDict()
        for i, residue_grp in enumerate(homology_info):

            if residue_grp["residue_chains"]:

                mcsa_id = residue_grp["mcsa_id"]
                if mcsa_id not in ordered_homology_info:
                    ordered_homology_info[mcsa_id] = {}

                reference_resid = residue_grp["residue_chains"][0]
                if reference_resid["is_reference"]:
                    ref_residue_id = reference_resid["resid"]
                    ref_pdb_id = reference_resid["pdb_id"]
                    ref_chain_name = reference_resid["chain_name"]
                    protein_key = ref_pdb_id + ref_chain_name
                else:
                    raise ValueError(
                        f"First residue was not the reference in mcsa_id {mcsa_id} residue number {i} !"
                    )

                if protein_key not in ordered_homology_info[mcsa_id]:
                    ordered_homology_info[mcsa_id][protein_key] = {}

                if ref_residue_id not in ordered_homology_info[mcsa_id][protein_key]:
                    ordered_homology_info[mcsa_id][protein_key][
                        ref_residue_id
                    ] = residue_grp

                # Now I want to figure out why!
                else:
                    uniprots_equal = (
                        residue_grp["residue_sequences"][0]["uniprot_id"]
                        and residue_grp["residue_sequences"][0]["uniprot_id"]
                        == ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                            "residue_sequences"
                        ][0]["uniprot_id"]
                    )
                    chains_equal = (
                        residue_grp["residue_chains"][0]["chain_name"]
                        == ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                            "residue_chains"
                        ][0]["chain_name"]
                    )
                    # is it because the same residue happens to be catalytic in different proteins?
                    if uniprots_equal and chains_equal:
                        # Here most differences between the two residue groups seem to stem from distinct roles
                        # Often the same residue will interact both through its side chain and its main chain etc. ec.
                        # In these cases we prefer the residue group with the side chain interaction role
                        if (
                            ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                                "function_location_abv"
                            ]
                            and not residue_grp["function_location_abv"]
                        ):
                            ordered_homology_info[mcsa_id][protein_key][
                                ref_residue_id
                            ] = residue_grp
                        elif (
                            not ordered_homology_info[mcsa_id][protein_key][
                                ref_residue_id
                            ]["function_location_abv"]
                            and residue_grp["function_location_abv"]
                        ):
                            pass

                        # it might be due to role differences. We add together all roles
                        elif (
                            ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                                "roles_summary"
                            ]
                            != residue_grp["roles_summary"]
                        ):
                            current_roles = set(
                                role.strip()
                                for role in ordered_homology_info[mcsa_id][protein_key][
                                    ref_residue_id
                                ]["roles_summary"].split(",")
                            )
                            new_roles = set(
                                role.strip()
                                for role in residue_grp["roles_summary"].split(",")
                            )
                            new_summary = ",".join(new_roles & current_roles)
                            ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                                "roles_summary"
                            ] = new_summary

                            # we deduplicate later
                            ordered_homology_info[mcsa_id][protein_key][ref_residue_id][
                                "roles"
                            ].extend(residue_grp["roles"])
                        else:
                            # Three cases are left but I dont really care.
                            # doesnt look lite it hurts me either
                            pass
                    else:
                        raise ValueError(
                            "This error should never be raised. It indicates serious annoation issues. Contact Antonio if it does."
                        )

            else:
                problematic_entries.add(mcsa_id)

        print(f"problematic entries are {sorted(list(problematic_entries))}")

        # sort outer keys
        ordered_homology_info = OrderedDict(sorted(ordered_homology_info.items()))

        # sort inner keys
        for mcsa_id, dict_of_reference_pdbchains in ordered_homology_info.items():
            for ref_pdbchain, residue_dict in dict_of_reference_pdbchains.items():
                residue_dict = OrderedDict(sorted(residue_dict.items()))
                ordered_homology_info[mcsa_id][ref_pdbchain] = residue_dict

        ################################################################################

        # start a new dict with the final format I want
        mcsa_pdb_resids: OrderedDict = OrderedDict()
        for mcsa_id, dict_of_reference_pdbchains in ordered_homology_info.items():

            if mcsa_id not in mcsa_pdb_resids:
                mcsa_pdb_resids[mcsa_id] = OrderedDict()

            for _, residue_dict in dict_of_reference_pdbchains.items():

                for index, residue in residue_dict.items():

                    for i, hom_residue in enumerate(residue["residue_chains"]):

                        # skip structures which are not reference or template structures
                        if hom_residue["is_reference"]:
                            pass
                        elif (
                            hom_residue["pdb_id"] + hom_residue["chain_name"]
                            in template_pdbchains
                        ):
                            pass
                        elif (
                            hom_residue["pdb_id"] + hom_residue["assembly_chain_name"]
                            in template_pdbchains
                        ):
                            pass
                        elif (
                            hom_residue["pdb_id"] + hom_residue["chain_name"][0]
                            in template_pdbchains
                        ):
                            pass
                        elif (
                            hom_residue["pdb_id"]
                            + hom_residue["assembly_chain_name"][0]
                            in template_pdbchains
                        ):
                            pass
                        else:
                            continue

                        if i == 0 and hom_residue["is_reference"]:
                            reference_pdbchain = (
                                hom_residue["pdb_id"] + hom_residue["chain_name"]
                            )
                        elif i == 0 and not hom_residue["is_reference"]:
                            raise ValueError("reference not found in first position!")

                        # TODO ask antonio about this!
                        # Fixing auth_resid and resid assignments in case of None???
                        if (
                            hom_residue["auth_resid"] is None
                            and hom_residue["resid"] is not None
                        ):
                            hom_residue["auth_resid"] = hom_residue["resid"]
                        elif (
                            hom_residue["auth_resid"] is not None
                            and hom_residue["resid"] is None
                        ):
                            hom_residue["resid"] = hom_residue["auth_resid"]

                        pdb_entry = HomologousPDB(
                            is_reference=hom_residue["is_reference"],
                            mcsa_id=mcsa_id,
                            pdb_id=hom_residue["pdb_id"],
                            reference_pdbchain=reference_pdbchain,
                            chain_name=hom_residue["chain_name"],
                            assembly_chain_name=hom_residue["assembly_chain_name"],
                            assembly=hom_residue["assembly"],
                            residues={},
                        )

                        # hopefuly this is sufficiently unique
                        pdb_identifier = pdb_entry.pdb_id + pdb_entry.chain_name

                        if pdb_identifier not in mcsa_pdb_resids[mcsa_id]:
                            mcsa_pdb_resids[mcsa_id][pdb_identifier] = pdb_entry

                        # trying to rescue some of the assembly information??
                        if (
                            mcsa_pdb_resids[mcsa_id][pdb_identifier].assembly is None
                            and hom_residue["assembly"] is not None
                        ):
                            mcsa_pdb_resids[mcsa_id][pdb_identifier].assembly = (
                                hom_residue["assembly"]
                            )

                        if (
                            index
                            not in mcsa_pdb_resids[mcsa_id][pdb_identifier].residues
                        ):
                            if hom_residue["is_reference"]:
                                mcsa_pdb_resids[mcsa_id][
                                    pdb_identifier
                                ].is_reference = True
                                mcsa_pdb_resids[mcsa_id][pdb_identifier].residues[
                                    index
                                ] = ReferenceCatalyticResidue(
                                    code=hom_residue["code"],
                                    resid=hom_residue["resid"],
                                    auth_resid=hom_residue["auth_resid"],
                                    function_location_abv=residue[
                                        "function_location_abv"
                                    ],
                                    ptm=residue["ptm"],
                                    roles_summary=[
                                        s.strip()
                                        for s in residue["roles_summary"].split(",")
                                    ],
                                    roles=list(
                                        role_dict["emo"]
                                        for role_dict in residue["roles"]
                                    ),
                                )
                            else:
                                mcsa_pdb_resids[mcsa_id][pdb_identifier].residues[
                                    index
                                ] = NonReferenceCatalyticResidue(
                                    code=hom_residue["code"],
                                    resid=hom_residue["resid"],
                                    auth_resid=hom_residue["auth_resid"],
                                    reference=(mcsa_id, reference_pdbchain, index),
                                )
                        else:
                            print(
                                f"the same index {index} already exists for the pdb {pdb_identifier} in mcsa {mcsa_id}"
                            )
                            # raise ValueError(f"the same index {index} already exists for the pdb {pdb_identifier} in mcsa {mcsa_id}")

        #################### Custom manual fixes #######################################
        # wrong auth_resid assignments in 1ddxA
        mcsa_pdb_resids[37]["1ddxA"].residues[172].auth_resid = 203  # GLN 172
        mcsa_pdb_resids[37]["1ddxA"].residues[176].auth_resid = 207  # HIS 176
        mcsa_pdb_resids[37]["1ddxA"].residues[353].auth_resid = 384  # LEU 353
        mcsa_pdb_resids[37]["1ddxA"].residues[354].auth_resid = 385  # TYR 354
        mcsa_pdb_resids[37]["1ddxA"].residues[357].auth_resid = 388  # HIS 357
        mcsa_pdb_resids[37]["1ddxA"].residues[495].auth_resid = 526  # GLY 495
        mcsa_pdb_resids[37]["1ddxA"].residues[499].auth_resid = 530  # SER 499

        # wrong auth_resid assignments in 3fr7A
        mcsa_pdb_resids[890]["3fr7A"].residues[181].auth_resid = 252  # LYS 181
        mcsa_pdb_resids[890]["3fr7A"].residues[244].auth_resid = 315  # ASP 244
        mcsa_pdb_resids[890]["3fr7A"].residues[248].auth_resid = 319  # GLU 248

        # missing residue 280 in 1esdA
        # is also not included in the table in the M-CSA!
        mcsa_pdb_resids[556]["1esdA"].residues[280] = NonReferenceCatalyticResidue(
            code="Trp", resid=280, auth_resid=280, reference=(556, "1escA", 280)
        )

        # templates in truth use the assembly_chain_id but too frequently they use
        # odd assembly chain names like AA instead of B etc. etc.
        # or often if assembly-chain-names are A-1, A-2 etc. templates have A, AA, ...
        # so I decided to use the unique chain_id instead
        # therefore we have to apply some fixes. I checked all of them manually
        mcsa_pdb_resids[11]["1qumD"] = mcsa_pdb_resids[11]["1qumA"]
        mcsa_pdb_resids[85]["1vmdAA"] = mcsa_pdb_resids[85]["1vmdB"]
        mcsa_pdb_resids[162]["1vasC"] = mcsa_pdb_resids[162]["1vasA"]
        mcsa_pdb_resids[182]["1oqfAA"] = mcsa_pdb_resids[182]["1oqfB"]
        mcsa_pdb_resids[375]["1otxCA"] = mcsa_pdb_resids[375]["1otxB"]
        mcsa_pdb_resids[404]["1asyC"] = mcsa_pdb_resids[404]["1asyA"]
        mcsa_pdb_resids[495]["1dmuB"] = mcsa_pdb_resids[495]["1dmuA"]
        mcsa_pdb_resids[759]["5ijgB"] = mcsa_pdb_resids[759]["5ijgA"]
        mcsa_pdb_resids[897]["1pviC"] = mcsa_pdb_resids[897]["1pviA"]

        # fix annotations in entry 638:
        metal_role = {
            "group_function": "metal ligand",
            "function_type": "interaction",
            "function": "metal ligand",
            "emo": "EMO_00116",
        }
        metal_summary = ["metal ligand"]
        mcsa_pdb_resids[638]["1i1iP"].residues[475].roles = (
            mcsa_pdb_resids[638]["1i1iP"].residues[503].roles
        )
        mcsa_pdb_resids[638]["1i1iP"].residues[475].roles_summary = (
            mcsa_pdb_resids[638]["1i1iP"].residues[503].roles_summary
        )

        mcsa_pdb_resids[638]["1i1iP"].residues[474].roles.append(metal_role)
        mcsa_pdb_resids[638]["1i1iP"].residues[503].roles = [metal_role]
        mcsa_pdb_resids[638]["1i1iP"].residues[478].roles.append(metal_role)
        mcsa_pdb_resids[638]["1i1iP"].residues[474].roles_summary = metal_summary
        mcsa_pdb_resids[638]["1i1iP"].residues[503].roles_summary = metal_summary
        mcsa_pdb_resids[638]["1i1iP"].residues[478].roles_summary = metal_summary

        mcsa_pdb_resids[585]["1lbuA"].residues[154].roles.append(metal_role)
        mcsa_pdb_resids[585]["1lbuA"].residues[161].roles.append(metal_role)
        mcsa_pdb_resids[585]["1lbuA"].residues[197].roles.append(metal_role)
        mcsa_pdb_resids[585]["1lbuA"].residues[154].roles_summary = metal_summary
        mcsa_pdb_resids[585]["1lbuA"].residues[161].roles_summary = metal_summary
        mcsa_pdb_resids[585]["1lbuA"].residues[197].roles_summary = metal_summary
        ################################################################################

        # NOTE:
        # Usually for a monometric biological assembly,
        # all the reference catalytic residues will be in one chain.
        # In this case 'is_reference' == True means
        # all the residues will be instances of 'ReferenceCatalyticResidue'
        # While 'is_reference' == False means
        # all residues will be instances of 'NonReferenceCatalyticResidue'

        # However, for mutlimeric catalytic sites,
        # the reference residues will also be split across two different chains.
        # In this case a HomologousPDB with 'is_reference' == True
        # may contain both reference and non-refererence residues.
        # If ANY of the residues in the chain are a 'ReferenceCatalyticResidue'
        # then the HomoogousPDB will be set with 'is_reference' == True

        ################################################################################

        # save the dict to a serializable json file
        serializable = {
            mcsa_id: {
                pdb_identifier: asdict(pdb_entry)
                for pdb_identifier, pdb_entry in pdb_entries.items()
            }
            for mcsa_id, pdb_entries in mcsa_pdb_resids.items()
        }
        with open(file, "w") as f:
            json.dump(serializable, f, indent=4, sort_keys=True)


def download_static_data(
    outdir: Path = Path("/home/ray/Documents/template_matching/enzymm/data"),
):
    get_mcsa_ec(outdir)
    get_mcsa_cath(outdir)
    check_other_files(outdir)
    make_pdb_sifts_df(outdir)
    download_mcsa_homolog_info(outdir)

    print("All files for annotation of jess results are there!")


if __name__ == "__main__":
    download_static_data()
