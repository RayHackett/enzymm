import io
import unittest

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

import pyjess

import enzymm
from . import test_data
from enzymm import template, jess_run
from enzymm.output import Tables, FullMatchTable, SimpleMatchTable, MatchResidueTable


class TestOutput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        template_files = [
            "jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/5_residues/results/csa3d_0045/csa3d_0045.cluster_1_1_1.2cxg_A227-A229-A257-A327-A328.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_2_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0421/csa3d_0421.cluster_1_1_2.1bf2_A229-A232-A230-A259-A375-A435-A510-A128.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0896/csa3d_0896.cluster_2_1_2.2qy1_A135-A179-A137-A230-A175.template.pdb",
        ]

        template_texts = []
        for path in template_files:
            with open(
                resource_files(enzymm).joinpath(path), "r"
            ) as f:  # ty:ignore[no-matching-overload]
                template_texts.append(f.read())

        template_list = []
        for txt in template_texts:
            template_list.append(template.AnnotatedTemplate.loads(txt))

        matcher = jess_run.Matcher(
            templates=template_list,
            filter_matches=True,
            skip_smaller_hits=False,
        )

        with resource_files(test_data).joinpath("1AMY.pdb").open() as f:
            molecule = pyjess.Molecule.load(f)  # ty:ignore[invalid-argument-type]

        cls.matches = matcher.run_single(molecule)

    def test_FullMatchTable(self):
        tbl = Tables.create(kind="full", matches=self.matches)

        self.assertIsInstance(tbl, FullMatchTable)
        self.assertEqual(tbl.columns(), FullMatchTable.columns())

        df = tbl.to_polars()
        # TODO write tests

        buffer = io.StringIO()
        tbl.write_tsv(file=buffer, header=False)
        with resource_files(test_data).joinpath("1amy_full_results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_SimpleMatchTable(self):
        tbl = Tables.create(kind="simple", matches=self.matches)

        self.assertIsInstance(tbl, SimpleMatchTable)
        self.assertEqual(tbl.columns(), SimpleMatchTable.columns())

        df = tbl.to_polars()
        # TODO write tests

        buffer = io.StringIO()
        tbl.write_tsv(file=buffer, header=False)
        with resource_files(test_data).joinpath("1amy_simple_results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_MatchResidueTable(self):
        tbl = Tables.create(kind="residue", matches=self.matches)

        self.assertIsInstance(tbl, MatchResidueTable)
        self.assertEqual(tbl.columns(), MatchResidueTable.columns())

        df = tbl.to_polars()
        # TODO write tests

        buffer = io.StringIO()
        tbl.write_tsv(file=buffer, header=False)
        with resource_files(test_data).joinpath("1amy_residue_results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())
