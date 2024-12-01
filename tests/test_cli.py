import unittest
import matcher
from pathlib import Path
import errno
import importlib.resources
from importlib.resources import files
from . import test_data
import io
from matcher._cli import main


class Test_CLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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
            "temp_result.tsv",
        ]
        cls.arguments_list = [
            "-l",
            list_path,
            "-t",
            selected_template_dir,
            "-o",
            "temp_result.tsv",
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
            "temp_result.tsv",
        ]

        cls.bad_argument_1 = [
            "-l",
            molecule_path,
            "-o",
            "temp_result.tsv",
        ]  # not a file in list file
        cls.bad_argument_2 = [
            "-l",
            list_path2,
            "-o",
            "temp_result.tsv",
        ]  # passing a dir in list file

        # TODO either remove the temp_result.tsv file or avoid writing it? I cant see how that oculd work though

    def test_default_main(self):
        self.assertEqual(main(self.arguments_normal), 0)
        self.assertEqual(main(self.arguments_list), 0)
        self.assertEqual(main(self.arguments_both), 0)

        retcode = main(self.bad_argument_1, stderr=io.StringIO())
        self.assertEqual(retcode, errno.ENOENT)
        retcode2 = main(self.bad_argument_2, stderr=io.StringIO())
        self.assertEqual(retcode2, errno.EISDIR)
