import numpy as np

def minimal_cover(elements_set, subsets, heuristic='default',
                  selected=(), extended_elements_set=None, depth=0):
    """Specific version of the "minimal subset cover greedy solution".

    This version allows for primary/extended coverage subsets (more doc later)

    Parameters
    ----------
    elements_set
      The set of all ements to cover

    subsets
      A list of (name, subset)

    heuristic
      A function ``((name, subset), selected) => value`` where ``name`` is the
      name of a subset, ``subset`` is what remains of the subset at this stage,
      ``selected`` is a list of already-selected subset names.

    selected
      (Recursion parameter, do not use.) Already-selected elements

    depth
      (Recursion parameter, do not use.). Depth of the recursion

    Returns
    -------

      None if no solution was found, else a collection of [(name, subset)...]
      in the order in which the subsets
    """
    if depth == 0:
        full_set = set().union(*[
            subset['extended']
            for name, subset in subsets
        ])
        if full_set != elements_set:
            raise ValueError('No full coverage solution exists !')
    if extended_elements_set is None:
        extended_elements_set = elements_set.union({})
    if (len(extended_elements_set) == 0) or (len(elements_set) == 0):
        return []
    subsets = list(subsets)
    subsets = [(n, s) for (n, s) in subsets if len(s['primary'])]

    name, subset = None, None
    for e in elements_set:
        containing_e = [(n, s) for (n, s) in subsets if e in s['primary']]
        if len(containing_e) == 1:
            name, subset = containing_e[0]
            ordered_subsets = subsets
            break

    if name is None:
        def sorting_heuristic(named_subset):
            if (heuristic == 'default'):
                return len(named_subset[1]['primary'])
            else:
                return heuristic(named_subset, selected)
        ordered_subsets = sorted(subsets, key=sorting_heuristic)

        name, subset = ordered_subsets.pop()
    primary, extended = subset['primary'], subset['extended']
    new_elements_set = elements_set.difference(primary)
    new_extended_elements_set = extended_elements_set.difference(extended)
    new_subsets = [
        (name_, {'primary': sub['primary'].difference(primary),
                 'extended': sub['extended'].difference(extended)})
        for (name_, sub) in ordered_subsets
    ]
    return [name] + minimal_cover(
        new_elements_set, new_subsets,
        heuristic=heuristic,
        selected=list(selected) + [subset],
        extended_elements_set=new_extended_elements_set,
        depth=depth + 1
)



def segments_to_array(segments, array_length):
    array = np.zeros(array_length)
    for start, end in segments:
        array[start: end] = 1
    return array

def group_overlapping_segments(segments, min_distance=10):
    if segments == []:
        return []
    returned_segments = [list(segments[0])]
    for start, end in segments[1:]:
        if start < returned_segments[-1][-1] + min_distance:
            if end > returned_segments[-1][-1]:
                returned_segments[-1][-1] = end
        else:
            returned_segments.append([start, end])
    return [tuple(s) for s in returned_segments]
