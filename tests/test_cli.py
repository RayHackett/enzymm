import unittest
import errno
import io
import tempfile
import shutil
from pathlib import Path

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

import enzymm
from enzymm._cli import main
from . import test_data


class Test_CLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tempfile = tempfile.NamedTemporaryFile(suffix=".tsv")
        molecule_path = str(
            Path(
                resource_files(test_data).joinpath(
                    "1AMY.pdb"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )
        list_path = str(
            Path(
                resource_files(test_data).joinpath(
                    "input_list.txt"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )
        list_path2 = str(
            Path(
                resource_files(test_data).joinpath(
                    "input_dir.txt"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )
        list_path3 = str(
            Path(
                resource_files(test_data).joinpath(
                    "input_list_bad.txt"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )
        list_path4 = str(
            Path(
                resource_files(test_data).joinpath(
                    "input_list_bad2.txt"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )

        selected_template_dir = str(
            Path(
                resource_files(enzymm).joinpath(
                    "jess_templates_20230210/5_residues/results/csa3d_0285/"
                )  # ty:ignore[invalid-argument-type]
            ).resolve()
        )

        cls.arguments_normal = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
        ]

        cls.arguments_simple = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "--simple-results",
        ]

        cls.arguments_simple_with_residues = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "--simple-results",
            "--per-residue-results",
        ]

        cls.arguments_simple_with_residues_parquet = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "--simple-results",
            "--per-residue-results",
            "--write-parquet",
        ]

        cls.arguments_normal_with_pdb = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "--pdbs",
            str(Path(tempfile.gettempdir(), "pdbs").resolve()),
            "--include-template",
        ]

        cls.arguments_normal_with_pdb_transformed = [
            "-i",
            molecule_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "--pdbs",
            str(Path(tempfile.gettempdir(), "pdbs_t").resolve()),
            "--include-template",
            "--transform",
        ]

        cls.arguments_unfiltered = [
            "-i",
            molecule_path,
            "-o",
            cls.tempfile.name,
            "-t",
            selected_template_dir,
            "-p",
            "2",
            "0.5",
            "0.5",
            "--unfiltered",
        ]

        cls.arguments_list = [
            "-l",
            list_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "-p",
            "5",
            "1.1",
            "1.5",
            "-c",
            "2",
        ]

        cls.arguments_both = [
            "-i",
            molecule_path,
            "-l",
            list_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
        ]

        cls.bad_argument_1 = [
            "-l",
            molecule_path,
            "-o",
            cls.tempfile.name,
        ]  # not a file in list file
        cls.bad_argument_2 = [
            "-l",
            list_path2,
            "-o",
            cls.tempfile.name,
        ]  # passing a dir in list file
        cls.bad_argument_3 = [
            "-o",
            cls.tempfile.name,
        ]  # no input

        cls.bad_argument_4 = [
            "-i",
            list_path,
            "-o",
            cls.tempfile.name,
        ]  # passing a something other than a PDB through -i

        cls.bad_argument_5 = [
            "-l",
            list_path3,
            "-o",
            cls.tempfile.name,
        ]  # passing a missing PDB in the list

        cls.bad_argument_6 = [
            "-l",
            list_path4,
            "-o",
            cls.tempfile.name,
        ]  # passing something which is not a PDB in the list

        cls.arguments_bad_unfiltered = [
            "-i",
            molecule_path,
            "-o",
            cls.tempfile.name,
            "-p",
            "2",
            "0.5",
            "0.5",
        ]  # expecting filtering with 0.5A

        cls.maxDiff = None

    @classmethod
    def tearDownClass(cls):
        cls.tempfile.close()
        shutil.rmtree(Path(tempfile.gettempdir(), "pdbs"))
        shutil.rmtree(Path(tempfile.gettempdir(), "pdbs_t"))

    def test_default_main(self):
        self.assertEqual(main(self.arguments_normal, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_unfiltered, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_list, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_both, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_simple, stderr=io.StringIO()), 0)
        self.assertEqual(
            main(self.arguments_simple_with_residues, stderr=io.StringIO()), 0
        )
        self.assertEqual(
            main(self.arguments_simple_with_residues_parquet, stderr=io.StringIO()), 0
        )

        retcode = main(self.bad_argument_1, stderr=io.StringIO())
        self.assertEqual(retcode, errno.ENOENT)
        retcode2 = main(self.bad_argument_2, stderr=io.StringIO())
        self.assertEqual(retcode2, errno.EISDIR)
        with self.assertRaises(ValueError):
            main(self.bad_argument_3, stderr=io.StringIO())
        with self.assertRaises(ValueError):
            main(self.bad_argument_4, stderr=io.StringIO())
        retcode3 = main(self.bad_argument_5, stderr=io.StringIO())
        self.assertEqual(retcode3, errno.ENOENT)
        with self.assertRaises(ValueError):
            main(self.bad_argument_6, stderr=io.StringIO())
        with self.assertRaises(ValueError):
            main(self.arguments_bad_unfiltered, stderr=io.StringIO())

    def test_writing_pdbs_query_ref(self):
        main(self.arguments_normal_with_pdb, stderr=io.StringIO())

        self.assertTrue(Path(tempfile.gettempdir(), "pdbs").is_dir())
        self.assertTrue(
            Path(Path(tempfile.gettempdir(), "pdbs"), "1AMY_matches.pdb").is_file()
        )

        with open(
            Path(Path(tempfile.gettempdir(), "pdbs"), "1AMY_matches.pdb"), "r"
        ) as f:
            actual = f.read()

        with resource_files(test_data).joinpath(
            "1AMY_matches_with_templates.pdb"
        ).open() as f:
            expected = f.read()

        self.assertMultiLineEqual(actual.strip(), expected.strip())

    def test_writing_pdbs_template_ref(self):
        main(self.arguments_normal_with_pdb_transformed, stderr=io.StringIO())

        self.assertTrue(Path(tempfile.gettempdir(), "pdbs_t").is_dir())
        self.assertTrue(
            Path(Path(tempfile.gettempdir(), "pdbs_t"), "1uh3_matches.pdb").is_file()
        )

        with open(
            Path(Path(tempfile.gettempdir(), "pdbs_t"), "1uh3_matches.pdb"), "r"
        ) as f:
            actual = f.read()

        with resource_files(test_data).joinpath(
            "1AMY_matches_with_templates_t.pdb"
        ).open() as f:
            expected = f.read()

        self.assertMultiLineEqual(actual.strip(), expected.strip())
