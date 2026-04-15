import enzymm
import unittest
import json
import tempfile
from pathlib import Path

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from enzymm import utils
from enzymm.jess_run import load_molecules
from . import test_data

try:
    from isal import igzip as gzip  # ty:ignore[unresolved-import]
except ImportError:
    import gzip  # ty:ignore[invalid-assignment]

try:
    import lz4.frame  # ty:ignore[unresolved-import]
except ImportError as err:
    lz4 = err  # ty:ignore[invalid-assignment]

try:
    import bz2
except ImportError as err:
    bz2 = err  # ty:ignore[invalid-assignment]

try:
    import lzma
except ImportError as err:
    lzma = err  # ty:ignore[invalid-assignment]


class TestUtils(unittest.TestCase):

    def test_ranked_argsort(self):
        self.assertEqual(utils.ranked_argsort([0, 4, 8, 6]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([-3, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 20, 9]), [1, 2, 4, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 20, 3, 20, 9]), [1, 4, 2, 4, 3])

    def test_chunks(self):
        self.assertEqual(list(utils.chunks(["a", "b", "c"], 2)), [("a", "b"), ("c",)])
        self.assertEqual(
            list(utils.chunks(["a", "b", "c"], 1)), [("a",), ("b",), ("c",)]
        )
        self.assertEqual(len(list(utils.chunks(set(["a", "b", "b", "c"]), 2))), 2)
        self.assertEqual(
            list(utils.chunks({1: "a", 2: "b", 3: "c"}, 2)), [(1, 2), (3,)]
        )

    # def test_request_url(self): no idea how to test this tbh

    def test_json_extract(self):
        self.assertEqual(
            utils.json_extract({1: "a", 2: "b", 3: {1: "c", 2: "d"}}, 2), ["b", "d"]
        )
        self.assertEqual(
            utils.json_extract([{1: "a", 2: "b"}, {1: "c", 2: "d"}], 2), ["b", "d"]
        )


class TestDecompression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.molecule_path = Path(
            resource_files(test_data).joinpath(
                "1AMY.pdb"
            )  # ty:ignore[invalid-argument-type]
        ).resolve()
        cls.molecule = load_molecules(molecule_paths=[cls.molecule_path])[0]

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_gzip_roundtrip(self):
        with open(self.molecule_path, "r", encoding="utf-8") as f:
            original = f.read()

        # using isal, write the original text to a temp compressed binary .gz file
        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY.pdb.gz"
        with gzip.open(path_to_tmp_file, "wb") as gz_out:
            gz_out.write(original.encode("utf-8"))

        loaded = load_molecules(molecule_paths=[path_to_tmp_file])[0]
        self.assertEqual(self.molecule, loaded)

    def test_bzip2_roundtrip(self):
        with open(self.molecule_path, "r", encoding="utf-8") as f:
            original = f.read()

        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY.pdb.bz2"
        with bz2.open(path_to_tmp_file, "wb") as bz2_out:
            bz2_out.write(original.encode("utf-8"))

        loaded = load_molecules(molecule_paths=[path_to_tmp_file])[0]
        self.assertEqual(self.molecule, loaded)

    def test_lz4_roundtrip(self):
        with open(self.molecule_path, "r", encoding="utf-8") as f:
            original = f.read()

        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY.pdb.lz4"
        with lz4.frame.open(path_to_tmp_file, "wb") as lz4_out:
            lz4_out.write(original.encode("utf-8"))

        loaded = load_molecules(molecule_paths=[path_to_tmp_file])[0]
        self.assertEqual(self.molecule, loaded)

    def test_lzma_roundtrip(self):
        with open(self.molecule_path, "r", encoding="utf-8") as f:
            original = f.read()

        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY.pdb.xz"
        with lzma.open(path_to_tmp_file, "wb") as xz_out:
            xz_out.write(original.encode("utf-8"))

        loaded = load_molecules(molecule_paths=[path_to_tmp_file])[0]
        self.assertEqual(self.molecule, loaded)

    def test_invalid_extension_raises(self):
        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY.pdb.zst"
        path_to_tmp_file.write_bytes(b"fake compressed data")

        with self.assertRaises((ValueError, OSError)):
            load_molecules(molecule_paths=[path_to_tmp_file])

    def test_corrupted_gz_raises(self):
        path_to_tmp_file = Path(self.tmpdir.name) / "1AMY_corrupt.pdb.gz"
        path_to_tmp_file.write_bytes(b"\x1f\x8b\x08\x00corrupted data here!!")

        try:
            from isal import igzip_lib

            expected_exc = igzip_lib.IsalError
        except ImportError:
            from gzip import BadGzipFile

            expected_exc = BadGzipFile

        with self.assertRaises(expected_exc):  # ty:ignore[invalid-argument-type]
            load_molecules(molecule_paths=[path_to_tmp_file])


class TestSetEncoder(unittest.TestCase):
    def test_encode_set(self):
        result = json.dumps({"numbers": {1, 2, 3, 4, 5}}, cls=utils.SetEncoder)
        self.assertEqual(result, '{"numbers": [1, 2, 3, 4, 5]}')
