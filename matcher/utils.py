import requests  # type: ignore
from requests.adapters import HTTPAdapter  # type: ignore
from urllib3.util.retry import Retry  # type: ignore
import typing
from typing import List, Tuple, Iterable, Iterator, TypeVar, Any, Literal
from itertools import islice
import time
import json


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


# this makes any set serializable. This allows me to write to json
# consider: https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
class SetEncoder(json.JSONEncoder):
    """
    Class to transform sets into lists to enable json serialization
    """

    def default(self, obj):
        if isinstance(obj, set):  # add more object types with elif isinstance etc.
            return list(obj)
        return super().default(obj)  # return parent class


# for running either single or multithread
class DummyPool:
    "Class to mimic a Threadpool. Used for single threaded runs."

    def map(self, function, iterable):
        return list(map(function, iterable))

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


T = TypeVar("T")


@typing.overload
def chunks(iterable: Iterable[T], n: Literal[2]) -> Iterator[Tuple[T, T]]: ...


@typing.overload
def chunks(iterable: Iterable[T], n: Literal[3]) -> Iterator[Tuple[T, T, T]]: ...


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
