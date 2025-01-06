import unittest
import matcher
from pathlib import Path
import errno
import importlib.resources
from importlib.resources import files
from . import test_data
import io
import tempfile
from matcher._cli import main


class Test_CLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tempfile = tempfile.NamedTemporaryFile(suffix=".tsv")
        molecule_path = str(Path(files(test_data).joinpath("1AMY.pdb")).resolve())

        list_path = str(Path(files(test_data).joinpath("input_list.txt")).resolve())
        list_path2 = str(Path(files(test_data).joinpath("input_dir.txt")).resolve())

        selected_template_dir = str(
            Path(
                importlib.resources.files(matcher).joinpath(
                    "jess_templates_20230210/5_residues/results/csa3d_0285/"
                )
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
        cls.arguments_list = [
            "-l",
            list_path,
            "-t",
            selected_template_dir,
            "-o",
            cls.tempfile.name,
            "-j",
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

    @classmethod
    def tearDownClass(cls):
        cls.tempfile.close()

    def test_default_main(self):
        self.assertEqual(main(self.arguments_normal, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_list, stderr=io.StringIO()), 0)
        self.assertEqual(main(self.arguments_both, stderr=io.StringIO()), 0)

        retcode = main(self.bad_argument_1, stderr=io.StringIO())
        self.assertEqual(retcode, errno.ENOENT)
        retcode2 = main(self.bad_argument_2, stderr=io.StringIO())
        self.assertEqual(retcode2, errno.EISDIR)
