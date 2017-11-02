"""All things sequencing simulation."""

import re
from collections import OrderedDict
from .SequencingMatches import (SequenceMatch, SequencingMatches,
                               SequencingMatchesSet)
from .biotools import reverse_complement

def simulate_sequencing(sequence, primers, read_length=600, gap=50,
                        linear=True, trim_nonmatching_primers=True):
    """

    Parameters
    ----------

    sequence

    primers

    read_length

    gap

    linear

    trim_nonmatching_primers

    Returns
    --------

    If primers is simply a Primer object, a SequencingMatches object
    featuring the primer matches on the sequence (`match.primer_matches`)
    and the associated read(s) (`match.read_matches`)

    If `primers` is an list/tuple of primers, an ordered dictionnary of the
    form {primer_name: SequencingMatches()}
    """

    if isinstance(primers, (list, tuple)):
        result = OrderedDict([
            (p.name, simulate_sequencing(sequence=sequence, primers=p,
                                         read_length=read_length,
                                         gap=gap, linear=linear))
            for p in primers
        ])
        if trim_nonmatching_primers:
            for name, matches in list(result.items()):
                if matches.read_matches == []:
                    result.pop(name)
        if len(result) == 0:
            return None
        return SequencingMatchesSet(result, linear=linear)

    primer = primers

    len_seq = len(sequence)
    record = sequence
    if not linear:
        sequence = sequence + sequence

    if hasattr(sequence, "seq"):
        sequence = str(record.seq)

    def simulate_one_read(reverse):

        seq = reverse_complement(sequence) if reverse else sequence

        primer_matches = [
            SequenceMatch(start=m.start(), end=m.end())
            for m in re.finditer(primer.sequence, seq)
        ]

        read_matches = [
            SequenceMatch(start=m.start + gap, end=m.start + gap + read_length)
            for m in primer_matches
        ]

        for match_group in (primer_matches, read_matches):
            for match in match_group:
                start, end = match.start, match.end
                if reverse:
                    start, end = len(seq) - end, len(seq) - start
                    match.strand = -1
                if start > len_seq:
                    if linear:
                        match_group.remove(match)
                    else:
                        start, end = start - len_seq, end - len_seq
                elif end < 0:
                    if linear:
                        match_group.remove(match)
                    else:
                        start, end = start + len_seq, end + len_seq
                match.start, match.end = start, end

        return primer_matches, read_matches

    primer_matches, read_matches = simulate_one_read(reverse=False)
    primer_matches_r, read_matches_r = simulate_one_read(reverse=True)
    all_primer_matches = primer_matches + primer_matches_r
    all_read_matches = read_matches + read_matches_r
    return SequencingMatches(record, all_primer_matches, all_read_matches)
