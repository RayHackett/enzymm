import requests # type: ignore
from requests.adapters import HTTPAdapter # type: ignore
from urllib3.util.retry import Retry # type: ignore
from pathlib import Path
import typing
from typing import Dict, Set, List, Tuple, Iterable, Iterator, TypeVar, Union, Any
from itertools import islice
import time
import json

def ranked_argsort(lst: List[int]) -> List[int]:
    """ Return a list of the same order in which the elements values correspond to their ranked values"""
    unique_values = sorted(set(lst))
    ranks = {v: i+1 for i, v in enumerate(unique_values)}
    return [ranks[i] for i in lst]

def request_url(url: str, acceptable_stati: List[int], timeout: int =10, max_retries: int =10):
    # TODO figure out maxretries

    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)

    r = session.get(url)
    request_counter = 0
    while request_counter < max_retries and r.status_code in range(500, 600): # server errors 500-599; 429 too many requests
        time.sleep(timeout)
        r = session.get(url)
        request_counter += 1
        
    if r.status_code == 429:
        time.sleep(int(r.headers["Retry-After"]))
        
    if r.status_code in acceptable_stati:
        if r.status_code == 200:
            return r
        else:
            return ''  
    else:
        print('Failed to get url', url, ' with status code: ', r.status_code)

# this makes any set serializable. This allows me to write to json
# consider: https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set): # add more object types with elif isinstance etc.
            return list(obj)
        return super().default(obj) # return parent class

X = TypeVar('X', str, int, float, bool)  # Add more types if needed
def json_extract(obj: Any, key: X) -> List[X]:
    """Recursively fetch values from nested dictionary"""
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

T = TypeVar('T')
def chunks(iterable: Iterable[T], n: int) -> Iterator[List[T]]:
    """Yield successive n-sized chunks from iterable."""
    iterable = iter(iterable)
    while chunk := list(islice(iterable, n)):
        yield chunk
