from typing import List, Tuple
import itertools
from operator import attrgetter
import math
import pyjess # type: ignore

__all__ = [
    "filter_hits"
]

def filter_hits(hits: List[pyjess.Hit], all_instances=False) -> List[pyjess.Hit]:
    # Given a list of hits, group that list into sub-lists per target template pair
    if not hits:
        raise ValueError('Did not recieve any hits')
    # For each sub-list find the best hit and append it to a list (of length equal to the number of sub-lists)
    # Group hit objects by the pair (hit.template, hit.target)
    def get_key(obj: pyjess.Hit) -> Tuple[str, str]:
        return obj.template.id, obj.molecule.id

    grouped_hits = [list(g) for _, g in itertools.groupby(sorted(hits, key=get_key), get_key)]

    best_hits = []
    for hits in grouped_hits:
        best_hits.extend(_best_match(hits))

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
        if not math.isclose(hit.determinant, 1, rel_tol=1e-05):
            continue # skips the rest of the loop
        # This is to discard self matches when cross-comparing templates
        # This bit is a bit messy because it actually uses the id of both objects so it could break if the user isnt careful
        # if hit.molecule.id != None:
        #     if hit.template.id.lower() == hit.molecule.id.lower():
        #         continue # skips the rest of the loop
        unique_hits.append(hit)

    if all_instances:  
        # do not filter by highest score
        return unique_hits # ! This is a list of hit objects
    else:
        if unique_hits:
            return [_get_highest_scored_hit(unique_hits)] # ! This returns a single hit so we turn it into a list
        else:
            return [] # return an empty list if all were self matches or had bad determinants
        