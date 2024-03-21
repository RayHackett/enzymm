from typing import List, Tuple
import itertools
from operator import attrgetter
import pyjess

__all__ = [
    "filter_hits"
]

def filter_hits(hits: List[pyjess.Hit], all_instances=False) -> List[pyjess.Hit]:
    # Given a list of hits, group that list into sub-lists per target template pair
    if not hits:
        raise ValueError('Did not recieve any hits')
    # For each sub-list find the best hit and append it to a list (of length equal to the number of sub-lists)
    # Group hit objects by the pair (hit.template, hit.target)
    def get_key(obj) -> Tuple[str, str]:
        return obj.template, obj.target

    grouped_hits = [list(g) for _, g in itertools.groupby(sorted(hits, key=get_key), get_key)]

    best_hits = []
    for hits in grouped_hits:
        best_hits.append(_best_match(hits))

    return best_hits

def _get_highest_scored_hit(hits: List[pyjess.Hit]) -> pyjess.Hit:
    # takes a list of matches and returns the match with the best (aka lowest) e-value score
    if not hits:
        raise ValueError('Did not recieve any unique hits')
    best_hit = min(hits, key=attrgetter('log_evalue'))
    return best_hit

def _best_match(hits: List[pyjess.Hit], all_instances=False) -> List[pyjess.Hit]:
    if not hits:
        raise ValueError('Did not recieve any hits')
    # takes a list of matches, filters out self-matches and returns the best match or matches
    unique_hits = []
    for hit in hits:
        # implement check for D-value == 1, otherwise continue
        if hit.determinant != 1.0:
            continue # skips the rest of the loop
        # This is to discard self matches when cross-comparing templates
        if hit.template == hit.target: # change this accoridng to above comments
            continue # skips the rest of the loop
        unique_hits.append(hit)
    if all_instances:  
        # do not filter by highest score
        return unique_hits # ! This is a list of hit objects
    else:
        return [_get_highest_scored_hit(unique_hits)] # ! This returns a single hit so we turn it into a list # TODO this is pretty bad handling honestly. I dont want a list in this case at all...
        