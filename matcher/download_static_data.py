from pathlib import Path
import json
import pandas as pd  # type: ignore
import urllib.request
import gzip
import shutil

from utils import request_url, SetEncoder  # type: ignore


def check_other_files():
    Path("./data").mkdir(parents=True, exist_ok=True)

    EC_cofactor_mapping = Path("./data/EClist_cofactors_forRH.csv")
    # AlphaFold_CATH = Path('./data/cath-v4_3_0.alphafold-v2.2022-11-22.tsv')
    PDBchain_to_CATH_Uniprot = Path("./data/pdb_chain_cath_uniprot.csv")
    PDBchain_to_EC_Uniprot = Path("./data/pdb_chain_enzyme.csv")
    PDBchain_to_InterPro = Path("./data/pdb_chain_interpro.csv")
    CATH_names_mapping = Path("./data/cath-domain-list.txt")

    if not Path(EC_cofactor_mapping).is_file():
        raise FileNotFoundError("EC to cofactor mapping csv from Neera is missing!")

    if not Path(PDBchain_to_CATH_Uniprot).is_file():
        urllib.request.urlretrieve(
            "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_cath_uniprot.csv.gz",
            "./data/pdb_chain_cath_uniprot.csv.gz",
        )

        with gzip.open("./data/pdb_chain_cath_uniprot.csv.gz", "rb") as f_gz:
            with open("./data/pdb_chain_cath_uniprot.csv", "wb") as f:
                shutil.copyfileobj(f_gz, f)
        # Mapping PDBchains to CATH name and Uniprot! Obtain via:
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

        Path.unlink("./data/pdb_chain_cath_uniprot.csv.gz")

    if not Path(PDBchain_to_EC_Uniprot).is_file():

        urllib.request.urlretrieve(
            "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_enzyme.csv.gz",
            "./data/pdb_chain_enzyme.csv.gz",
        )

        with gzip.open("./data/pdb_chain_enzyme.csv.gz", "rb") as f_gz:
            with open("./data/pdb_chain_enzyme.csv", "wb") as f:
                shutil.copyfileobj(f_gz, f)

        Path.unlink("./data/pdb_chain_enzyme.csv.gz")

        # Mapping PDBchains to EC and Uniprot
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

    if not Path(PDBchain_to_InterPro).is_file():

        urllib.request.urlretrieve(
            "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_interpro.csv.gz",
            "./data/pdb_chain_interpro.csv.gz",
        )

        with gzip.open("./data/pdb_chain_interpro.csv.gz", "rb") as f_gz:
            with open("./data/pdb_chain_interpro.csv", "wb") as f:
                shutil.copyfileobj(f_gz, f)

        Path.unlink("./data/pdb_chain_interpro.csv.gz")

        # Mapping PDBchains to EC and Uniprot
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

    if not Path(CATH_names_mapping).is_file():
        # Mapping CATH ids to CATH numbers and CATH domain names
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt",
            "./data/cath-domain-list.txt",
        )

    # if not Path(AlphaFold_CATH).is_file():
    #     # Mapping UniProt AFDB Proteins to CATH Domains
    #     urllib.request.urlretrieve('ftp://orengoftp.biochem.ucl.ac.uk/alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv', './data/cath-v4_3_0.alphafold-v2.2022-11-22.tsv')


def make_pdb_sifts_df():
    if not Path("./data/pdb_sifts.json").is_file():
        check_other_files()

        PDBchain_to_CATH_Uniprot = Path("./data/pdb_chain_cath_uniprot.csv")
        PDBchain_to_EC_Uniprot = Path("./data/pdb_chain_enzyme.csv")
        PDBchain_to_InterPro = Path("./data/pdb_chain_interpro.csv")
        CATH_names_mapping = Path("./data/cath-domain-list.txt")

        # load them all as df
        pdb_cath_uniprot_df = pd.read_csv(
            PDBchain_to_CATH_Uniprot, comment="#", dtype=str, na_filter=False
        )
        pdb_enzyme_uniprot_df = pd.read_csv(
            PDBchain_to_EC_Uniprot, comment="#", dtype=str, na_filter=False
        )
        pdb_interpro_uniprot_df = pd.read_csv(
            PDBchain_to_InterPro, comment="#", dtype=str, na_filter=False
        )
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
        pdb_interpro_uniprot_df["PDBCHAIN"] = (
            pdb_interpro_uniprot_df["PDB"] + pdb_interpro_uniprot_df["CHAIN"]
        )

        sifts_dict = {}
        for row in pdb_enzyme_uniprot_df.itertuples():
            if row.PDBCHAIN not in sifts_dict:
                sifts_dict[row.PDBCHAIN] = {
                    "uniprot_id": row.ACCESSION,
                    "cath": set(),
                    "ec": set(),
                    "interpro": set(),
                }
            if row.EC_NUMBER != "?":
                sifts_dict[row.PDBCHAIN]["ec"].add(row.EC_NUMBER)
        for row in pdb_cath_uniprot_df.itertuples():
            if row.PDBCHAIN not in sifts_dict:
                sifts_dict[row.PDBCHAIN] = {
                    "uniprot_id": row.SP_PRIMARY,
                    "cath": set(),
                    "ec": set(),
                    "interpro": set(),
                }
            sifts_dict[row.PDBCHAIN]["cath"].add(row.CATH_NUMBER)
        for row in pdb_interpro_uniprot_df.itertuples():
            if row.PDBCHAIN not in sifts_dict:
                sifts_dict[row.PDBCHAIN] = {
                    "uniprot_id": None,
                    "cath": set(),
                    "ec": set(),
                    "interpro": set(),
                }
            sifts_dict[row.PDBCHAIN]["interpro"].add(row.INTERPRO_ID)

        # merge the dataframes
        Path("./data").mkdir(parents=True, exist_ok=True)
        with open(Path("./data/pdb_sifts.json"), "w") as f:
            json.dump(sifts_dict, f, indent=4, sort_keys=True, cls=SetEncoder)


def get_mcsa_ec():
    # m-csa entry map to EC number
    # every m-csa entry only has one EC number
    if not Path("./data/MCSA_EC_mapping.json").is_file():
        ec_dict = {}
        for page in range(1, 12):
            url = (
                "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page="
                + str(page)
            )
            general = request_url(url, [200])
            general_page = general.json()

            for i in range(len(general_page["results"])):
                mcsa_id = int(general_page["results"][i]["mcsa_id"])
                if mcsa_id not in ec_dict:
                    ec_dict[mcsa_id] = ""
                if len(general_page["results"][i]["all_ecs"]) > 1:
                    print(mcsa_id)
                for ec_num in general_page["results"][i]["all_ecs"]:
                    ec_dict[mcsa_id] = ec_num

            Path("./data").mkdir(parents=True, exist_ok=True)
            with open(Path("./data/MCSA_EC_mapping.json"), "w") as f:
                json.dump(ec_dict, f, indent=4, sort_keys=True)


def get_mcsa_cath():
    # m-csa entry map to CATH id
    # entry have multiple CATH ids
    if not Path("./data/MCSA_CATH_mapping.json").is_file():
        mcsa_cath = {}
        for page in range(1, 12):
            url = (
                "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page="
                + str(page)
            )
            general = request_url(url, [200])
            general_page = general.json()

            for i in range(len(general_page["results"])):
                cath_ids = set()
                for resid in general_page["results"][i]["residues"]:
                    for chain in resid["residue_chains"]:
                        cath_ids.add(chain["domain_cath_id"])

                if "" in cath_ids:
                    cath_ids.remove("")
                cath_ids = list(cath_ids)
                mcsa_cath[general_page["results"][i]["mcsa_id"]] = cath_ids

            Path("./data").mkdir(parents=True, exist_ok=True)
            with open(Path("./data/MCSA_CATH_mapping.json"), "w") as f:
                json.dump(mcsa_cath, f, indent=4, sort_keys=True)


def download_static_data():
    get_mcsa_ec()
    get_mcsa_cath()
    check_other_files()
    make_pdb_sifts_df()

    print("All files for annotation of jess results are there!")


if __name__ == "__main__":
    download_static_data()
