"""Class for one-line sequencing report generation from a sequencing batch"""

import os
from collections import defaultdict
import flametree

import matplotlib.pyplot as plt
from Bio import SeqIO
from flametree import file_tree
from .ReadReferenceMatches import SequencingRead, ReadReferenceMatchesSet
from .Primer import Primer


class Source:

    def __call__(self, name):
        return self.get(self.sanitize_name(name), None)

    @staticmethod
    def sanitize_name(name):
        if "." in name:
            name = name.split(".")[0]
        return name.lower().replace(" ", "_").replace("-", "_")


class PrimersFastaSource(Source):
    """Primer source using a fasta file for primers names and sequences."""

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.primers_dict = {
            self.sanitize_name(primer.name): primer
            for primer in Primer.list_from_fasta(fasta_file)
        }
        self.get = self.primers_dict.get

class PrimersSpreadsheetSource(Source):
    """Primer source using a spreadsheet for primers names and sequences."""

    def __init__(self, spreadsheet):
        self.spreadsheet = spreadsheet
        if isinstance(spreadsheet, str):
            primers_list = Primer.list_from_spreadsheet(filepath=spreadsheet)
        else:
            primers_list = Primer.list_from_spreadsheet(dataframe=spreadsheet)
        self.primers_dict = {
            self.sanitize_name(primer.name): primer
            for primer in primers_list
        }
        self.get = self.primers_dict.get

class ConstructsFolderSource(Source):
    """Loads all files in the folders as Biopython records.
    """

    def __init__(self, dirs, extensions=("gb", "gbk"), file_format="genbank"):
        self.records = {}
        for c_folder in dirs:
            root = flametree.file_tree(c_folder)
            self.records.update({
                self.sanitize_name(f._name_no_extension):
                    SeqIO.read(f._path, file_format)
                for f in root._all_files
                if f._extension in extensions
            })
        print (len(self.records))
        self.get = self.records.get


class SequencingReportGenerator:
    """Reads a series of sequencing reads from a zip. Make a report"""

    def __init__(self, primers_source=None, constructs_source=None,
                 default_linearity=False):
        self.primers_source = primers_source
        self.constructs_source = constructs_source
        self.default_linearity = default_linearity

    def get_read_infos(self, filename):
        """Return clone_id, primer_name, construct_name from the filename."""
        return filename.split("_")[:2]

    def classify_reads(self, reads):
        """Return a dict {construct_id: [associated_reads]}."""
        classified_reads = defaultdict(lambda *a: [])
        for read in reads:
            clone_id, primer_name, construct_name = \
                self.get_read_infos(read.read_name)
            read.primer = self.primers_source(primer_name)

            classified_reads[clone_id].append(read)
        for clone, seqlist in classified_reads.items():
            seqlist.sort(key=lambda r: r.primer.name)
        self.classified_reads = classified_reads
        return classified_reads

    def plot_matches_set(self, matches_set, filepath, title=None):
        ax = matches_set.plot(
            plot_coverage=True,
            plot_reference=True, reference_ax=None,
            figsize="auto", features_filters=(),
            features_properties=None, reference_reads_shares="auto")
        if title is not None:
            ax.set_title(title)
        ax.figure.savefig(filepath, format="png", bbox_inches="tight")
        plt.close(ax.figure)

    def make_report(self, ab1_files=None, ab1_zip_file=None,
                    target="my_folder", perc_identity=100,
                    plot_params=None, replace=True):
        """"""
        if ab1_zip_file is not None:
            ab1_files = [
                f for f in file_tree(ab1_files)._all_files
                if f._extension == "ab1"
            ]
        reads = [SequencingRead.from_ab1_file(f) for f in ab1_files]
        reads = [r for r in reads if set(r.read_sequence) != set('N')]

        classified_reads = self.classify_reads(reads)
        errors = []
        records = []

        root = file_tree(target, replace=replace)

        for clone_id, reads in classified_reads.items():
            _, _, construct_name = self.get_read_infos(reads[0].read_name)
            construct = self.constructs_source(construct_name)
            if construct is None:
                errors.append(construct_name + ": unknown construct")
                continue
            linear = construct.__dict__.get('linear', self.default_linearity)
            matches_set = ReadReferenceMatchesSet.from_reads(
                construct, reads, perc_identity=perc_identity, linear=linear
            )
            self.plot_matches_set(matches_set, title=clone_id,
                                  filepath=root._file(clone_id + ".png"))
            matches_set.to_genbank(root._file(clone_id + ".gb"))
            records += [
                dict(
                    read=read_name,
                    clone_id=clone_id,
                    construct=construct_name,
                    primer=matches.primer.name,
                    primer_position=matches.primer_start,
                    farthest_reading_span=matches.farthest_reading_span,
                    longest_match=matches.longest_match_size,
                    total_matches_length=matches.total_matches_length,
                    average_read_quality=matches.read.average_quality,
                    read_length=len(matches.read.read_qualities)
                )
                for read_name, matches in
                sorted(matches_set.read_reference_matches.items())
            ]
        columns = ["read", "construct", "primer", "primer_position",
                   "longest_match", "total_matches_length",
                   "average_read_quality", "read_length"]
        csv_content = "\n".join(
            [",".join(columns)] + [
                ",".join([('' if (record[col] is None) else str(record[col]))
                          for col in columns])
                for record in records
            ]
        )
        root._file("report.csv").write(csv_content)
        return root._close(), errors
