"""Bla."""
import tempfile
import subprocess
import time
import os
import sys
import zipfile
from copy import deepcopy

import numpy as np
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

PYTHON3 = (sys.version_info[0] == 3)

if PYTHON3:
    from io import StringIO, BytesIO
    StringBytesIO = BytesIO

else:
    from StringIO import StringIO
    BytesIO = StringIO
    StringBytesIO = StringIO

def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.
    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.
    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def blast_sequences(sequences=None, fasta_file=None,
                    blast_db=None, subject=None, word_size=4,
                    perc_identity=80, num_alignments=1000, num_threads=3,
                    use_megablast=True, evalue=None, ungapped=True):
    """Return a Biopython BLAST record of the given sequence BLASTed
    against the provided database.

    Parameters
    ----------

    sequences
      Either an ATGC string or a list of ATGC strings or a dict {name: seq:}

    subject
      Either a path to a fasta (.fa) file or an ATGC string. Subject to blast
      against.

    word_size
      Word size to use in the blast

    perc_identity
      Minimal percentage of identical nucleotides in a match for it to be kept

    num_alignments
      Number of alignments to keep

    num_threads
      Number of threads for the BLAST

    use_megablast
      Whether to use Megablast.

    ungapped
      No-gaps matches only ?



    Examples
    --------

    >>> blast_record = blast_sequence("ATTGTGCGTGTGTGCGT", "blastdb/ecoli")
    >>> for alignment in blast_record.alignments:
    >>>     for hit in alignment.hsps:
    >>>         print (hit.identities)
    """

    if isinstance(sequences, str):
        sequences = [sequences]

    if isinstance(sequences, (list, tuple)):
        sequences = {"seq_%d" % i: seq for i, seq in enumerate(sequences)}

    xml_file, xml_name = tempfile.mkstemp(".xml")
    fasta_file, fasta_name = tempfile.mkstemp(".fa")
    with open(fasta_name, "w+") as f:
        for (name, seq) in sequences.items():
            f.write(">%s\n%s\n" % (name, seq))

    remove_subject = True
    close_subject = False
    if subject is not None:
        close_subject = True
        if isinstance(subject, str):
            if subject.endswith(".fa"):
                remove_subject = False
            else:
                subject = [subject]
        if isinstance(subject, (list, tuple)):
            subject = {"subject_%d" % i: seq for i, seq in enumerate(subject)}
        if isinstance(subject, dict):
            subject_file, fasta_subject_name = tempfile.mkstemp(".fa")
            with open(fasta_subject_name, "w+") as f:
                for (name, seq) in subject.items():
                    f.write(">%s\n%s\n" % (name, seq))
            subject = fasta_subject_name
    else:
        close_subject = False
    p = subprocess.Popen([
        "blastn", "-out", xml_name,
        "-outfmt", "5",
        "-num_alignments", str(num_alignments),
        "-query", fasta_name] +
        (["-db", blast_db] if blast_db is not None
         else ['-subject', subject]) +
        (["-ungapped"] if ungapped else []) +
        (["-evalue", str(evalue)] if evalue else []) +
        (["-task", "megablast"] if use_megablast else []) + [
        "-word_size", str(word_size),
        "-num_threads", str(num_threads),
        "-perc_identity", str(perc_identity),
        "-dust", "no"
    ], close_fds=True)
    res, blast_err = p.communicate()
    error = None
    for i in range(3):
        try:
            with open(xml_name, "r") as f:
                res = list(NCBIXML.parse(f))
        except ValueError as err:
            error = err
            time.sleep(0.1)
        else:
            break
    else:
        raise ValueError("Problem reading the blast record: " + str(error))
    for j in range(3):
        try:
            os.fdopen(xml_file, 'w').close()
            os.fdopen(fasta_file, 'w').close()
            os.remove(xml_name)
            os.remove(fasta_name)
            if close_subject:
                open(subject, 'w').close()
                if remove_subject:
                    os.remove(subject)
        except IOError as err:
            error = err
            time.sleep(0.1)
        else:
            break
    return res


# RETIRING THIS ONE ASAP:
def read_records_from_zip(zip_path):
    """Return SeqRecords from all FASTA/GENBANK files in the zip."""
    with zipfile.ZipFile(zip_path, 'r') as archive:
        extensions_types = {".ab1": "abi", ".abi": "abi", ".gb": "genbank",
                            ".gbk": "genbank", ".fa": "fasta",
                            ".fasta": "fasta"}
        extract = {}
        failed_files = []
        for f in archive.filelist:
            name, ext = os.path.splitext(f.filename)
            try:
                if ext in extensions_types:
                    content = StringBytesIO(archive.read(f.filename))
                    extract[f.filename] = SeqIO.read(content,
                                                     extensions_types[ext])
            except:
                failed_files.append(f.filename)
    return extract, failed_files

def rotate_circular_record(record, n_bases):
    """Changes the starting point of a circular SeqRecord by n_bases bases."""
    new_record = deepcopy(record)
    new_record.seq = record.seq[n_bases:] + record.seq[:n_bases]
    for f in new_record.features:
        f.location += (-n_bases)
        if max(f.location.start, f.location.end) <= 0:
            f.location += len(record)
    return new_record

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



def get_segment_coordinates(center, segment_length, sequence_length):
    """Return max(0, c - s/2) - min(L, c + L/2).

    Where c=center, s=segment_length, L=sequence_length.
    """
    half = int(segment_length / 2)
    start = max(0, min(center - half,  sequence_length - segment_length))
    end = start + segment_length
    return start, end

def find_best_primer_locations(sequence, size_range=(15, 25),
                               tm_range=(55, 70)):
    """Quickly compute all overhangs in the sequence.

    This function uses the heuristic {A, T}=2degC, {G, C}=4degC to compute
    melting temperatures.

    This function uses vectorial operations for speed. The results are also
    cached.
    """
    lmin, lmax = size_range
    tmin, tmax = tm_range

    table = np.zeros((lmax + 1 - lmin, len(sequence)))
    cumsum = np.cumsum([4 if nuc in "GC" else 2 for nuc in sequence])
    for i, oh_size in enumerate(range(lmin, lmax + 1)):
        arr = cumsum[oh_size:] - cumsum[:-oh_size]
        start = int(oh_size / 2)
        end = start + len(arr)
        table[i, start:end] = arr
        table[i, :start] = table[i, start]
        table[i, end:] = table[i, end-1]
    scores = - (table - tmin) * (table - tmax)
    best_sizes_indices = scores.argmax(axis=0)
    best_sizes = lmin + best_sizes_indices
    validities = np.choose(best_sizes_indices, scores) >= 0
    osizes_and_validities = zip(best_sizes, validities)
    return [
      None if not valid
      else get_segment_coordinates(i, ovh_size, len(sequence))
      for i, (ovh_size, valid) in enumerate(osizes_and_validities)
    ]

def find_non_unique_segments(sequence, perc_identity=80):
    blast_record = blast_sequences(sequence, subject=sequence,
                                   perc_identity=perc_identity,
                                   ungapped=False, word_size=4)[0]
    segments_with_alignments = sorted(set([
        (h.query_start, h.query_end)
        for al in blast_record.alignments
        for h in al.hsps
        if (h.query_start, h.query_end) != (1, len(sequence))
    ]))
    return group_overlapping_segments(segments_with_alignments)

def load_record(filename, linear=True, name='auto'):
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    record.linear = linear
    if name == 'auto':
        name = os.path.splitext(os.path.basename(filename))[0]
    record.id = name
    record.name = name.replace(" ", "_")[:20]
    return record

def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                    margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )
