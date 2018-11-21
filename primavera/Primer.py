from Bio import SeqIO
from .biotools import reverse_complement
import pandas

class Primer:
    """Class representing a sequencing primer

    Parameters
    ----------

    name
      Name (label) of the primer

    sequence
      An ATGC string

    metadata
      A dictionnary {field, value}
    """

    # TODO: try this, see if it makes things faster:
    # slots_= ['name', 'sequence', 'metadata']

    def __init__(self, name=None, sequence=None, metadata=None):
        """Initialize"""
        self.name = name
        self.sequence = sequence
        self.metadata = {} if metadata is None else metadata

    def find_location(self, sequence):
        """Return the (start, end, strand) of the primer in the sequence.

        The convention is that always (start < end), the strand indicates
        the direction of the homology.
        """
        ind = sequence.find(self.sequence)
        strand = 1
        if ind == -1:
            ind = sequence.find(reverse_complement(self.sequence))
            if ind == -1:
                return None
            else:
                ind = len(sequence) - ind - len(self.sequence)
                strand = -1
        start, end = ind, ind + len(self.sequence)
        return start, end, strand

    def __repr__(self):
        return "%s-%s" % (self.name, self.sequence)

    @staticmethod
    def list_from_fasta(fasta_file):
        """Return a list of Primer objects, by reading a FASTA file."""
        records = SeqIO.parse(fasta_file, 'fasta')
        return [
            Primer(name=rec.name, sequence=str(rec.seq))
            for rec in records
        ]

    @staticmethod
    def list_to_fasta(primers_list, filepath=None):
        """Write a FASTA file from a list of Primer objects."""
        fasta = "\n\n".join([
            ">%s\n%s" % (primer.name, primer.sequence)
            for primer in primers_list
        ])
        if filepath is not None:
            with open(filepath, 'w') as f:
                f.write(fasta)
        return fasta

    @staticmethod
    def list_from_spreadsheet(dataframe=None, filepath=None):
        """Return a list of Primer objects from a spreadsheet.

        The spreadsheet should have columns name/sequence."""
        if dataframe is None:
            if filepath.lower().endswith('csv'):
                dataframe = pandas.read_csv(filepath)
            else:
                dataframe = pandas.read_excel(filepath)

        return [
            Primer(name=r['name'], sequence=r['sequence'],
                   metadata=dict((k, v) for (k, v) in r.items()
                                 if k not in ['name', 'sequence']))
            for r in dataframe.to_dict(orient='record')
        ]

    @staticmethod
    def list_to_spreadsheet(primers_list, filepath=None):
        """Convert a list of Primer objects to a spreadsheet."""

        dataframe = pandas.DataFrame.from_records(
            columns=['name', 'sequence'] + sorted(set(
                field
                for primer in primers_list
                for field in primer.metadata
            )),
            data=[
                dict([('name', primer.name), ('sequence', primer.sequence)] +
                     list(primer.metadata.items()))
                for primer in primers_list
            ]
        )

        if filepath is not None:
            if filepath.lower().endswith('csv'):
                dataframe.write_csv(filepath)
            else:
                dataframe.write_excel(filepath)
        return dataframe
