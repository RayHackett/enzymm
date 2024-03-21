# common functions script

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pathlib import Path
from typing import List, Iterator
import time

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

def convert_sets_to_lists(obj):
    # This function recursively turns sets to lists inside nested dictionaries
    if isinstance(obj, set):
        return list(obj)
    elif isinstance(obj, dict):
        return {key: convert_sets_to_lists(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_sets_to_lists(item) for item in obj]
    else:
        return obj

def fasta_batching(input_file, splits, out_path, out_name):
    with open(input_file, 'r') as f:
        data = f.read()

    seq_list = data.split('\n>') # also removes the >

    new_seq_list = [seq_list[0]]
    for seq in seq_list[1:]:
        new_seq_list.append('>' + seq) # after the first we want to add > back

    # split Uniprot_list into batches of size splits
    batches = [new_seq_list[i:i + splits] for i in range(0, len(new_seq_list), splits)]       

    batch_list =[]
    for i, chunk in enumerate(batches):
        batch_number = i + 1
        batch_name = Path(out_path, out_name + str(batch_number) + '.fasta')
        with open(batch_name, 'w') as f:
            f.write('\n'.join(chunk))
        batch_list.append(batch_name)

    return batch_list

def file_batching(input_file, splits, out_path, out_name):
    with open(input_file, 'r') as f:
        data = f.read()
    data_list = data.split()

    # split Uniprot_list into batches of size splits
    batches = [data_list[i:i + splits] for i in range(0, len(data_list), splits)]

    batch_list =[]
    for i, chunk in enumerate(batches):
        batch_number = i + 1
        batch_name = Path(out_path, out_name + str(batch_number) + '.txt')
        with open(batch_name, 'w') as f:
            f.write('\n'.join(chunk))
        batch_list.append(batch_name)

    return batch_list

"""Extract nested values from a JSON tree."""
def json_extract(obj, key):
    """Recursively fetch values from nested JSON."""
    arr = []

    def extract(obj, arr, key):
        """Recursively search for values of key in JSON tree."""
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

def chunks(lst: List, n: int) -> Iterator[List]:
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

