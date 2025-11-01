from __future__ import annotations

from pathlib import Path
import json
import pandas as pd  # type: ignore
import csv
from typing import List, Tuple, Set, Dict, Optional, ClassVar, Any
from dataclasses import dataclass

from template_info import get_list_of_template_pdbchains

# This script requires files fron enzymm/data
# pdb_sifts.json
# catalytic_residue_homologs_information
# EClist_cofactors_forRH.csv

__all__ = [
    "Query",
]


@dataclass
class Query:
    """Class for storing Query annotation info derived from external data sources
    Required Annotations are uniprot, pdb_id, chain, cath ids, ec numbers, catalytic_site annotations from M-CSA
    Optional are Cofacor annotations, InterPro Annotations"""

    pdbchain: str
    pdb_id: str
    chain_id: str
    uniprot_id: str
    ec: List[str]
    cath: List[str]
    cofactors: Optional[List[str]]

    _PDB_SIFTS: ClassVar[Dict[str, Dict[str, List[str]]]]
    _homologs_info: ClassVar[Dict[str, Any]]
    _cofactor_dict: ClassVar[Dict[str, Set[str]]]

    @classmethod
    def load(cls, pdbchain: str) -> Query:
        """Function to create a Query object from a pdbchain string"""

        pdb_id = pdbchain[0:4].lower()
        chain_id = pdbchain[4:]  # case sensetive

        # superseded entries
        if pdb_id == "5esy":
            subdict = cls._PDB_SIFTS.get("8f9y" + chain_id)
        elif pdb_id == "7ayl":
            subdict = cls._PDB_SIFTS.get("7zsz" + chain_id)

        else:
            # read uniprot, ec and cath from sifts mappings
            subdict = cls._PDB_SIFTS.get(pdb_id + chain_id)

        if subdict is not None:
            uniprot = str(subdict.get("uniprot_id"))
            ec_list = list(subdict.get("ec"))  # type: ignore
            cath_list = list(subdict.get("cath"))  # type: ignore

        else:
            # TODO this is an ugly workaround:
            # Technically this is a bug either with missing annotations in SIFTS
            # or due to wierd chain name conventions in .cif files
            for chain_symbol in ["A", "B", "C"]:
                subdict = cls._PDB_SIFTS.get(pdb_id + chain_symbol)
                if subdict is None:
                    uniprot = ""
                    ec_list = []
                    cath_list = []

                else:
                    uniprot = str(subdict.get("uniprot_id"))
                    ec_list = list(subdict.get("ec"))  # type: ignore
                    cath_list = list(subdict.get("cath"))  # type: ignore
                    break

            if subdict is None:
                print(
                    f"Missing SIFTS annoations for {pdb_id+chain_id}. Likely there is no Uniprot ID"
                )

        # read in cofactors from Neeras list
        cofactors = set()
        for ec in ec_list:
            if ec in cls._cofactor_dict:
                cofactors.update(cls._cofactor_dict[ec])

        return Query(
            pdbchain=pdb_id + chain_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            uniprot_id=uniprot,
            ec=ec_list,
            cath=cath_list,
            cofactors=list(cofactors),
        )

    @property
    def mcsa_ids(self) -> List[int]:
        subdict = self._homologs_info[self.pdb_id]
        return list(subdict.get("mcsa_ids"))

    @property
    def auth_catalytic_sites(self) -> List[Tuple[int, str]]:

        subdict = self._homologs_info[self.pdb_id].get(self.chain_id)
        if subdict is None:
            subdict = self._homologs_info[self.pdb_id].get("A")
        if subdict is None:
            raise ValueError(f"Missing info for {self.pdb_id, self.chain_id}")

        return list(subdict.get("auth_catalytic_residues"))

    @property
    def assigned_catalytic_sites(self) -> List[Tuple[int, str]]:

        subdict = self._homologs_info[self.pdb_id].get(self.chain_id)
        if subdict is None:
            subdict = self._homologs_info[self.pdb_id].get("A")
        if subdict is None:
            raise ValueError(f"Missing info for {self.pdb_id, self.chain_id}")

        return list(subdict.get("assigned_catalytic_residues"))

    @property
    def is_reference(self) -> bool:
        return self._homologs_info[self.pdb_id]["is_reference"]

    # @cached_property
    # def interpro(self) -> Dict[str, List[str]]:
    #     url = 'https://www.ebi.ac.uk/interpro/api/entry/interpro/structure/pdb/' + self.pdb_id

    #     interpro_domains = set()
    #     interpro_families = set()
    #     interpro_superfam = set()
    #     r = request_url(url, [200, 204])
    #     if r:
    #         data = r.json()
    #         # by scraping only the stuff associated with certain amino acid positions this would be more accurate
    #         # look not only at the metadata stuff. but that gets a bit more complicated
    #         for i in data['results']:
    #             if i['metadata']['type'] == 'domain':
    #                 interpro_domains.add(i['metadata']['accession'])
    #             elif i['metadata']['type'] == 'family':
    #                 interpro_families.add(i['metadata']['accession'])
    #             elif i['metadata']['type'] == 'homologous_superfamily':
    #                 interpro_superfam.add(i['metadata']['accession'])
    #     else:
    #         print('No InterPro Response for:')
    #         print(url)

    #     return {'interpro_domains': list(interpro_domains), 'interpro_families': list(interpro_families), 'interpro_superfamilies': list(interpro_superfam)}

    # @cached_property
    # def catalytic(self):
    #     has_catalytic_go = False
    #     try:
    #         if self.uniprot_id != 'None':
    #             url = 'https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28accession%3A{}%29%20AND%20%28go%3A0003824%29%29'.format(self.uniprot_id)
    #             r = request_url(url, [200])
    #             if r.text:
    #                 has_catalytic_go = True
    #     except:
    #         pass

    #     return has_catalytic_go


############# Loading Files with annotation data ######################
outdir = Path("../../enzymm/data/")

# Source: CATH, EC and InterPro from PDB-SIFTS through mapping to the pdbchain
# This is now part of the template_matching tool, inside /matcher/data/
with open(Path(outdir, "pdb_sifts_new.json"), "r") as f:
    Query._PDB_SIFTS = json.load(f)

# Source: M-CSA downloads
with open(Path(outdir, "catalytic_residue_homologs_information.json"), "r") as f:
    homologs = json.load(f)
    Query._homologs_info = {}
    for mcsa_id, mcsa_dict in homologs.items():
        for pdbchain, pdbchain_dict in mcsa_dict.items():
            name = pdbchain_dict["pdb_id"].lower()

            chain = pdbchain_dict["chain_name"]

            if name not in Query._homologs_info:
                Query._homologs_info[name] = {
                    "is_reference": bool(pdbchain_dict["is_reference"]),
                    "mcsa_ids": set(),
                }
            Query._homologs_info[name]["mcsa_ids"].add(mcsa_id)
            if chain not in Query._homologs_info[name]:
                Query._homologs_info[name][chain] = {
                    "auth_catalytic_residues": set(),
                    "assigned_catalytic_residues": set(),
                }

            all_res_dict = pdbchain_dict["residues"]

            for res_idx, res_dict in all_res_dict.items():
                auth_resname = str(res_dict["auth_resid"])
                resname = str(res_dict["resid"])
                code = res_dict["code"].upper()

                if "None" not in auth_resname:
                    Query._homologs_info[name][chain]["auth_catalytic_residues"].add(
                        (int(auth_resname), code)
                    )
                if "None" not in resname:
                    Query._homologs_info[name][chain][
                        "assigned_catalytic_residues"
                    ].add((int(resname), code))

    # Properly assign assembly chain names to pdb chain names:
    Query._homologs_info["1qum"]["D"] = Query._homologs_info["1qum"]["A"]
    Query._homologs_info["1vmd"]["AA"] = Query._homologs_info["1vmd"]["B"]
    Query._homologs_info["1vas"]["C"] = Query._homologs_info["1vas"]["A"]
    Query._homologs_info["1oqf"]["AA"] = Query._homologs_info["1oqf"]["B"]
    Query._homologs_info["1otx"]["CA"] = Query._homologs_info["1otx"]["B"]
    Query._homologs_info["1asy"]["C"] = Query._homologs_info["1asy"]["A"]
    Query._homologs_info["1dmu"]["B"] = Query._homologs_info["1dmu"]["A"]
    Query._homologs_info["5ijg"]["B"] = Query._homologs_info["5ijg"]["A"]
    Query._homologs_info["1pvi"]["C"] = Query._homologs_info["1pvi"]["A"]
    Query._homologs_info["6l3p"]["BA"] = Query._homologs_info["6l3p"]["B"]
    Query._homologs_info["1vmd"]["BA"] = Query._homologs_info["1vmd"]["B"]
    Query._homologs_info["1rah"]["CB"] = Query._homologs_info["1rah"]["C"]

# I got this csv from Neera originally. Data is from 2020
# This is now part of the template_matching tool, inside /matcher/data/
with open(Path(outdir, "EClist_cofactors_forRH.csv"), "r") as f:
    cofactor_df = pd.read_csv(f)
    Query._cofactor_dict = {}
    for index, row in cofactor_df.iterrows():
        if row["EC"] not in Query._cofactor_dict:
            Query._cofactor_dict[row["EC"]] = set()
        Query._cofactor_dict[row["EC"]].add(row["Cof_ID"])


def write_query_tsv(outdir, query_list):
    # TODO fix writing for object attributes/properties that are lists
    with open(
        Path(outdir, "template_pdb_query_info.tsv"), "w", newline="", encoding="utf-8"
    ) as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "pdbchain",
                "pdb_id",
                "chain_id",
                "uniprot_id",
                "query_mcsa_id",
                "query_ec",
                "query_cath",
                "query_is_reference",
                "auth_catalytic_sites",
                "assigned_catalytic_sites",
                # "interpro_domains",
                # "interpro_families",
                # "interpro_superfamilies",
                # "catalytic",
                "cofactors",
            ]
        )
        for query in query_list:
            try:
                writer.writerow(
                    [
                        query.pdbchain,
                        query.pdb_id,
                        query.chain_id,
                        query.uniprot_id,
                        query.mcsa_ids,
                        query.ec,
                        query.cath,
                        query.is_reference,
                        query.auth_catalytic_sites,
                        query.assigned_catalytic_sites,
                        # query.interpro['interpro_domains'],
                        # query.interpro['interpro_families'],
                        # query.interpro['interpro_superfamilies'],
                        # query.catalytic,
                        query.cofactors,
                    ]
                )
            except KeyError as e:
                raise ValueError(
                    f"Issue with template pdbchain {query.pdbchain}"
                ) from e

    # TODO what should i retrun? a path to the output file? a bool if any hits were found?
    # return


if __name__ == "__main__":

    print(outdir.resolve())

    template_pdbchains_list = get_list_of_template_pdbchains(separator="")

    query_list = [Query.load(i) for i in template_pdbchains_list]

    write_query_tsv(outdir=outdir, query_list=query_list)
