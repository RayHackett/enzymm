from typing import (
    Dict,
    List,
    Tuple,
    Iterable,
    Iterator,
    TypeVar,
    Any,
    Callable,
    Literal,
    BinaryIO,
    TextIO,
    Union,
)
from itertools import islice
import json
import contextlib
import io
import pathlib

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

T = TypeVar("T")
U = TypeVar("U")


def ranked_argsort(lst: List[int]) -> List[int]:
    """
    Return a list of the same order in which the elements values correspond to their ranked values

    Arguments:
        lst: List on which to operate

    Returns:
        `List[int]`: List of value ranks of each element in order.
    """
    unique_values = sorted(set(lst))
    ranks = {v: i + 1 for i, v in enumerate(unique_values)}
    return [ranks[i] for i in lst]


# this makes any set serializable. This allows me to write to json
# consider: https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
class SetEncoder(json.JSONEncoder):
    """
    Class to transform sets into lists to enable json serialization
    """

    def default(self, obj):  # ty:ignore[invalid-method-override]
        if isinstance(obj, set):  # add more object types with elif isinstance etc.
            return list(obj)
        return super().default(obj)  # return parent class


# for running either single or multithread
class DummyPool:
    "Class to mimic a Threadpool. Used for single threaded runs."

    def map(self, function: Callable[[T], U], iterable: Iterable[T]) -> List[U]:
        return list(map(function, iterable))

    def imap(self, function: Callable[[T], U], iterable: Iterable[T]) -> Iterable[U]:
        return map(function, iterable)

    def imap_unordered(
        self, function: Callable[[T], U], iterable: Iterable[T]
    ) -> Iterable[U]:
        return map(function, iterable)

    def starmap(self, function, list_of_iterables):
        return list(function(*args) for args in list_of_iterables)

    def __enter__(self):
        return self

    def __exit__(self, exc_value, exc_type, traceback):
        return False


X = TypeVar("X", str, int, float, bool)  # Add more types if needed


def json_extract(obj: Any, key: X) -> List[X]:
    """
    Recursively fetch values from nested dictionary.

    Will also enter lists of dictionaries.

    Arguments:
        obj: Object instance from which to fetch values
        key: String, Int, float or bool key for which associated values will be fetched
    """
    arr: List[X] = []

    def extract(obj: Any, arr: List[X], key: X):
        """Recursively search for values of key in nested dictionary or list tree"""
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k == key:
                    arr.append(v)
                elif isinstance(v, (dict, list)):
                    extract(v, arr, key)
        elif isinstance(obj, list):
            for item in obj:
                extract(item, arr, key)
        return arr

    values = extract(obj, arr, key)
    return values


def chunks(iterable: Iterable[T], n: int) -> Iterator[Tuple[T, ...]]:
    """
    Yield successive n-sized chunks from iterable.

    Arguments:
        iterable: Iterable to chunk
        n: Number of elements in each chunk
    """
    iterable = iter(iterable)
    while chunk := tuple(islice(iterable, n)):
        yield chunk


PROTEINOGENIC_AMINO_ACIDS: Dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

TEMPLATE_ATOM_SELECTIOM: Dict[str, Tuple[str, str, str]] = {
    "ALA": ("C", "CA", "CB"),
    "ARG": ("CZ", "NH1", "NH2"),
    "ASN": ("CG", "OD1", "ND2"),
    "ASP": ("CG", "OD1", "OD2"),
    "CYS": ("CA", "CB", "SG"),
    "GLN": ("CD", "OE1", "NE2"),
    "GLU": ("CD", "OE1", "OE2"),
    "GLY": ("O", "C", "CA"),
    "HIS": ("CG", "CD2", "ND1"),
    "ILE": ("CA", "CB", "CG1"),
    "LEU": ("CA", "CB", "CG"),
    "LYS": ("CD", "CE", "NZ"),
    "MET": ("CG", "SD", "CE"),
    "PHE": ("CE1", "CZ", "CE2"),
    "PRO": ("O", "C", "CA"),
    "SER": ("CA", "CB", "OG"),
    "THR": ("CA", "CB", "OG1"),
    "TRP": ("NE1", "CZ2", "CH2"),
    "TYR": ("CE1", "CZ", "OH"),
    "VAL": ("CA", "CB", "CG1"),
    "ANY": ("O", "C", "CA"),
    "PTM": ("C", "CA", "CB"),
}

MATCH_MODE_TYPE = Literal[
    -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 100, 101, 102, 103, 104, 105, 106, 107
]

MATCH_MODE = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 100, 101, 102, 103, 104, 105, 106, 107]

PTMS_NAMED_IN_TEMPLATES = [
    "LLP",
    "TPQ",
    "SEC",
    "KCX",
    "NEP",
    "FGP",
    "SEP",
    "SMC",
    "CSS",
]

# some common special residues in pdb structures
SPECIAL_AMINO_ACIDS = [
    "ASX",
    "GLX",
    "LLP",
    "TPQ",
    "SEC",
    "PYL",
    "KCX",
    "NEP",
    "FGP",
    "UNK",
    "MSE",
    "SEP",
    "SMC",
    "TPO",
    "PTR",
    "HYP",
    "CME",
    "CSS" "CSO",
    "CSD",
    "PCA",
    "MLY",
    "DAL",
    "DAR",
    "DSG",
    "ORN",
    "PTM",
]

_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"


@contextlib.contextmanager
def zopen(path: Union[str, pathlib.Path, BinaryIO]) -> Iterator[TextIO]:
    """Open a file with optional compression in binary mode."""
    with contextlib.ExitStack() as ctx:
        if isinstance(path, (str, pathlib.Path)):
            file = ctx.enter_context(open(path, "rb"))
        else:
            file = ctx.enter_context(
                io.BufferedReader(path)  # ty:ignore[invalid-argument-type]
            )
        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            file = ctx.enter_context(gzip.open(file, mode="rt"))
        elif peek.startswith(_BZ2_MAGIC):
            if isinstance(bz2, ImportError):
                raise RuntimeError(
                    "File compression is LZMA but lzma is not available"
                ) from lz4
            file = ctx.enter_context(bz2.open(file, mode="rt"))
        elif peek.startswith(_XZ_MAGIC):
            if isinstance(lzma, ImportError):
                raise RuntimeError(
                    "File compression is LZMA but lzma is not available"
                ) from lz4
            file = ctx.enter_context(lzma.open(file, mode="rt"))
        elif peek.startswith(_LZ4_MAGIC):
            if isinstance(lz4, ImportError):
                raise RuntimeError(
                    "File compression is LZ4 but python-lz4 is not installed"
                ) from lz4
            file = ctx.enter_context(lz4.frame.open(file, mode="rt"))
        yield file
