"""Class to automatically select available and new primers for sequencing runs.
"""

from copy import deepcopy
import re
from collections import defaultdict

import pandas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

import flametree
from proglog import TqdmProgressBarLogger, ProgressBarLogger
from dnachisel import (AvoidPattern, RepeatedKmerPattern, HomopolymerPattern,
                       DnaOptimizationProblem)

from .sequencing_simulation import simulate_sequencing
from .biotools import (reverse_complement, find_non_unique_segments,
                       find_best_primer_locations)
from .tools import minimal_cover, segments_to_array, group_overlapping_segments
from .Primer import Primer


class PrimerSelectorLogger(TqdmProgressBarLogger):
    """Custom logger class adapted to the logger selector."""
    def __init__(self, bars=('record', 'primer'), notebook='default'):
        ignored_bars = set(('record', 'primer')).difference(bars)
        TqdmProgressBarLogger.__init__(self, bars=bars, notebook=notebook,
                                       ignored_bars=ignored_bars,
                                       min_time_interval=0.2)

class PrimerSelector:
    """A selector to compute the best primers to sequence a set of constructs.

    Examples
    --------

    >>> selector = PrimerSelector()
    >>> selected_primers = selector.select_primers(records, available_primers)
    >>> selector.plot_coverage(records, selected_primers, 'my_report.pdf'

    Parameters
    -----------

    read_range
      The experimentally measured range (start, end) so that, when a primer
      anneals in the sequence at index i, the range ``[i + start, i + end]`` will
      be correctly sequenced.

    size_range
      Size range (min, max) for the size of the primers, in nucleotides.

    tm_range
      Acceptable melting temperature range for the primers (in Celsius), as
      computed using the self.read_rangec A/T=2C, G/C=4C

    primer_conditions
      A list of functions of the form ``primer_sequence => True/False``.
      Primers for which at least one condition returns False will not be
      considered.

    primer_reuse_bonus
      Weight that the availability of the primer should have in the decision
      to select this primer. A higher value of this parameter leads to
      solutions where a higher less new primershave to be ordered, but more
      sequencing reactions have to be done. Set to e.g. 200 to test if there
      exists solutions involving solely already-available primers.

    logger
      Leave to 'bars' for default progress-bar logger, to None for no logger,
      or any Proglog ProgressBarLogger object.

    coverage_resolution
      When the user provides a record with "cover" features to indicate where
      to cover, the coverage points used by the algorithm are the 1-in-N
      nucleotides along the feature region, where N is this parameter.

    nucleotide_resolution
      If above 1, only every N primers will be considered when listing all
      the potential new primers (one around each nucleotide), where N is this
      number.
    """

    def __init__(self, read_range=(150, 800), size_range=(16, 25),
                 tm_range=(55, 70), primer_conditions=(),
                 primer_reuse_bonus=2, logger='bars',
                 homology_percentage=80,
                 nucleotide_resolution=1,
                 coverage_resolution=5):
        self.read_range = read_range
        self.size_range = size_range
        self.tm_range = tm_range
        self.homology_percentage = homology_percentage
        self.primers_conditions = primer_conditions
        self.nucleotide_resolution = nucleotide_resolution
        self.coverage_resolution = coverage_resolution
        self.primer_reuse_bonus = 2
        if logger == 'bars':
            logger = PrimerSelectorLogger()
        if logger is None:
            logger = ProgressBarLogger()
        self.logger = logger
        self.extension = coverage_resolution + 1

    def select_primers(self, records, available_primers=(),
                       strand='any', new_primers_prefix='P',
                       new_primers_digits=6):
        """Select primers to sequence the given records.

        Parameters
        ----------
        records
          A list of biopython records to sequence. The zones to cover in the
          record should be indicated by a feature of type ``misc_feature`` and
          label ``cover``. The zones where no primers are desired should be
          indicated by a feature of type ``misc_feature`` and label
          ``no_primer``.

        available_primers
          List of Primer objects representing the available primers.

        new_primers_prefix
          Prefix to use for names of the new primers

        new_primers_digits
          The new primers will have names of the form P000435, with a number
          of digits provided by this parameter.

        Returns
        -------
        selected_primers
          A list of lists of primers, one list of primers for each consrtruct.

        """
        if strand == 'both':
            reverse_records = []
            for record in records:
                new_record = record.reverse_complement()
                new_record.linear = record.linear
                new_record.id = record.id
                reverse_records.append(new_record)

            primers_sets = []
            for strand in self.logger.iter_bar(strand=['5-3', '3-5']):
                extra_prefix, recs = {'3-5': ('_r', reverse_records),
                                      '5-3': ('', records)}[strand]
                primers_sets.append(self.select_primers(
                    recs,
                    available_primers=available_primers,
                    strand='5-3',
                    new_primers_prefix=new_primers_prefix + extra_prefix,
                    new_primers_digits=new_primers_digits
                ))

            forward_primers, reverse_primers = primers_sets
            return [
                forward + rev
                for forward, rev in zip(forward_primers, reverse_primers)
            ]

        available_primers_dict = {p.sequence: p for p in available_primers}
        available_primers_seqs = set([p.sequence for p in available_primers])

        # COMPUTE PRIMERS AND COVERAGES
        indices_to_cover = {}
        primers_coverages = defaultdict(lambda *a: {'primary': set(),
                                                    'extended': set()})
        self.logger(message='Analyzing the records...')
        for record in self.logger.iter_bar(record=records):
            indices_to_cover[record.id] = {
                ind: '%s_%03d' % (record.id, ind)
                for ind in self.compute_indices_to_cover(record)
            }
            coverages = self.compute_all_primers_coverage_on_record(
                record, available_primers=available_primers_seqs,
                indices_to_cover=indices_to_cover[record.id],
                strand=strand)
            for primer, cov in coverages.items():
                primers_coverages[primer]['primary'].update(cov['primary'])
                primers_coverages[primer]['extended'].update(cov['extended'])
        # FIND GLOBAL MINIMAL COVER
        self.last_coverages = coverages
        self.logger(message='Selecting primers, please wait...')

        elements_set = set(
            index
            for rec_id, named_indices in indices_to_cover.items()
            for index in named_indices.values()
        )
        if len(elements_set) == 0:
            raise ValueError(
                "The provided problem has no specified locations "
                "to cover. Make sure you specify which parts of the sequence "
                "to cover, e.g. by adding annotations with label `cover` and "
                "type misc_feature in your records.")

        def heuristic(named_subset, selected):
            name, subset = named_subset
            primer_is_reused = name in available_primers_seqs
            min_index = min([
                int(re.match(r'.*_(\d*)', name_).groups()[0])
                for name_ in subset['primary']
            ])
            reuse_bonus = self.primer_reuse_bonus * primer_is_reused
            return (len(subset['primary']) + reuse_bonus, -min_index)

        subsets = primers_coverages.items()
        primers_cover = minimal_cover(elements_set, subsets=subsets,
                                      heuristic=heuristic)
        # print(list(primers_coverages.keys()))
        # REORGANIZE AND NAME THE SELECTED PRIMERS
        available_primers_names = [p.name for p in available_primers]
        selected_primers = []
        selected_primer_from_seq = {}
        for primer_seq in primers_cover:
            if primer_seq in available_primers_dict:
                name = available_primers_dict[primer_seq].name
                primer = available_primers_dict[primer_seq]
                infos = primer.metadata.get('infos', '')
                meta = {'available': True, 'infos': infos}
                meta.update(primer.metadata)
                primer = Primer(primer.name, primer.sequence, metadata=meta)
            else:
                name = self.generate_primer_name(
                    available_primers_names=available_primers_names,
                    prefix=new_primers_prefix,
                    n_digits=new_primers_digits
                )

                try:
                    part, index = self.find_subsequence_in_records(
                        sequence=primer_seq, records=records)
                    infos = "From part %s (at position %d)" % (part, index)
                except ValueError:
                    infos = "No containing part could be identified"
                available_primers_names.append(name)
                primer = Primer(name, primer_seq,
                                metadata={'available': False, 'infos': infos})
            selected_primers.append(primer)
            selected_primer_from_seq[primer_seq] = primer

        # CHOOSE A MINIMAL PRIMER COVER FOR EACH CONSTRUCT
        per_record_primers = []
        self.logger(message='Finding minimal covers for each construct...')
        for record in self.logger.iter_bar(record=records):
            elements = set(indices_to_cover[record.id].values())
            subcovers = {
                prim_seq: {
                    'primary': primers_coverages[prim_seq]['primary']
                               .intersection(elements),
                    'extended': primers_coverages[prim_seq]['extended']
                                .intersection(elements)
                }
                for prim_seq in primers_cover
            }
            subsets = deepcopy(list(subcovers.items()))
            # subsets = list(subcovers.items())
            sub_primers_cover = minimal_cover(elements, subsets=subsets,
                                              heuristic=heuristic)

            # All primers selected for this construct, sorted by availability.
            sub_selected_primers = sorted([
                selected_primer_from_seq[primer_seq]
                for primer_seq in sub_primers_cover
            ], key=lambda p: p.metadata['available'])

            per_record_primers.append(sub_selected_primers)

        return per_record_primers

    def compute_indices_to_cover(self, record):
        """List all indices in the record which should be covered.

        These indices are equidistant points inside the user-defined zones to
        cover in the record.

        The use determines the zones to cover via features of
        type ``misc_feature`` and label 'cover'.
        """
        segments_to_cover = [
            sorted([int(f.location.start), int(f.location.end)])
            for f in record.features
            if f.location is not None
            and f.type == 'misc_feature'
            and "".join(f.qualifiers.get('label', '')) == 'cover'
        ]
        res = self.coverage_resolution
        return set([
            int(indice) % len(record)
            for (start, end) in segments_to_cover
            for indice in np.linspace(start, end, 1.0 * (end - start) / res)
        ])

    @staticmethod
    def locate_primer_sequence(primer, sequence):
        """Find the location (start, end, strand) of a primer in the sequence.

        Return None if the primer sequence and its reverse complement are not
        found in the sequence.
        """
        ind = sequence.find(primer)
        strand = 1
        if ind == -1:
            ind = sequence.find(reverse_complement(primer))
            if ind == -1:
                return None
            else:
                strand = -1
        start, end = ind, ind + len(primer)
        return start, end, strand

    def compute_forbidden_patterns_locations(self, record):
        """Return an array where ``arr[i] == 1`` means that i is surrounded by
        a user-forbidden pattern."""
        pattern_constraints = [AvoidPattern(HomopolymerPattern(c, 5))
                               for c in 'ATGC']
        kmer_constraints = [
            AvoidPattern(RepeatedKmerPattern(n_repeats, kmer_size))
            for n_repeats, kmer_size in [(4, 2), (3, 3), (2, 4)]
        ]
        problem = DnaOptimizationProblem(
            sequence=record,
            constraints=pattern_constraints + kmer_constraints
        )
        constraints_breaches = group_overlapping_segments([
            (f.location.start, f.location.end)
            for ev in problem.constraints_evaluations()
            for f in ev.locations_to_biopython_features()
            if not ev.passes
        ])
        return segments_to_array(constraints_breaches, len(record))

    def compute_user_forbidden_locations(self, record):
        """Return an array where ``arr[i] == 1`` means that i is surrounded by
        a user-forbidden location."""
        forbidden_segments = [
            sorted([int(f.location.start), int(f.location.end)])
            for f in record.features
            if f.location is not None
            and f.type == 'misc_feature'
            and "".join(f.qualifiers.get('label', '')) == 'no_primer'
        ]
        return segments_to_array(forbidden_segments, len(record))

    def compute_nonunique_segments_locations(self, record):
        """Return an array where ``arr[i] == 1`` means that i is surrounded by
        a non-unique location."""
        sequence = str(record.seq)
        non_unique_segments = find_non_unique_segments(
            sequence, perc_identity=self.homology_percentage)
        return segments_to_array(non_unique_segments, len(record))

    def compute_all_forbidden_locations(self, record):
        """Return an array indicating which positions should be avoided.

        We take into account forbidden patterns, user-forbidden locations,
        and non-unique locations.
        ``arr[i] == 1`` indicates that position i should be avoided in the
        record.

        """
        return np.maximum(*(f(record) for f in (
            self.compute_nonunique_segments_locations,
            self.compute_forbidden_patterns_locations,
            self.compute_user_forbidden_locations
        )))

    def compute_sequence_primers(self, record, strand='any'):
        """Return primers for the sequence, one around each index.

        The primers are chosen to fit the length and melting temperature
        specified by the class parameters.

        Parameters
        ----------
        record
          The record in which to list the primers

        Returns
        -------
        primers
          List of primers sequences.

        """
        sequence = str(record.seq)
        L = len(sequence)
        rev_sequence = reverse_complement(sequence)
        locations = find_best_primer_locations(
            sequence,
            tm_range=self.tm_range,
            size_range=self.size_range
        )
        if self.nucleotide_resolution > 1:
            locations = locations[::self.nucleotide_resolution]
        result = []
        if strand in ['5-3', 'any']:
            result += [
                (sequence[l[0]: l[1]], (l[0], l[1], 1))
                for l in locations if l is not None
            ]
        if strand in ['3-5', 'any']:
            result += [
                (rev_sequence[L - l[1]: L - l[0]], (l[0], l[1], -1))
                for l in locations if l is not None
            ]
        return list(set(result))

    def compute_all_primers_coverage_on_record(self, record, indices_to_cover,
                                               available_primers,
                                               strand='any'):
        """Return, for each primer, the list of indices covered in this record.

        Parameters
        ----------
        record
          The record in which to search

        indices_to_cover
          List of indices to cover (defined by the user) in this record.

        available_primers
          List of candidate primers.

        Returns
        -------
        primers_coverage
          ``{primer_name: list_of_indices_covered_in_this_record}``.

        """
        sequence = str(record.seq)
        linear = record.__dict__.get('linear', True)
        L = len(sequence)
        forbidden_locations = self.compute_all_forbidden_locations(record)
        sequence_primers = self.compute_sequence_primers(record, strand=strand)
        reusable_primers_sequences = [
            (primer, location) for (primer, location) in [
                (primer, self.locate_primer_sequence(primer, sequence))
                for primer in available_primers
            ]
            if (location is not None) and
            ((strand == 'any') or
             ((location[2] == 1) and (strand == '5-3')) or
             ((location[2] == -1) and (strand == '3-5')))
        ]
        primers = reusable_primers_sequences + sequence_primers
        iter_bar = self.logger.iter_bar
        primers_coverages = {}
        for primer, (start, end, strand) in iter_bar(primer=primers):
            if forbidden_locations[start:end].max():
                continue
            if not all(cond(primer) for cond in self.primers_conditions):
                continue

            read_range_min = self.read_range[0] + self.extension
            read_range_max = self.read_range[1] - self.extension

            if strand == 1:
                coverage_start = end + read_range_min
                coverage_end = end + read_range_max
            else:
                coverage_start = start - read_range_max
                coverage_end = start - read_range_min
            if (not linear) and (coverage_start < 0):
                a, b, c, d = -1000, coverage_end, L + coverage_start, L + 1000
            elif (not linear) and (coverage_end > L):
                a, b, c, d = -1000, coverage_end - L, coverage_start, L + 1000
            else:
                a, b, c, d = (coverage_start, coverage_end,
                              coverage_start, coverage_end)
            a_ext = a - self.extension
            b_ext = b + self.extension
            c_ext = c - self.extension
            d_ext = d + self.extension
            primers_coverages[primer] = dict()
            primers_coverages[primer]['extended'] = set([
                ind
                for ind in indices_to_cover
                if (a_ext <= ind <= b_ext) or (c_ext <= ind <= d_ext)
            ])
            primers_coverages[primer]['primary'] = set([
                indices_to_cover[ind]
                for ind in primers_coverages[primer]['extended']
                if (a <= ind <= b) or (c <= ind <= d)
            ])
            primers_coverages[primer]['extended'] = set([
                indices_to_cover[ind]
                for ind in primers_coverages[primer]['extended']
            ])
            # print(len(primers_coverages[primer]['primary']))

        return primers_coverages

    def generate_primer_name(self, prefix='P', available_primers_names=(),
                             n_digits=6):
        """Return a suitable primer name, considering existing primer names.

        The result will be of the form P000425 where 'P' is the prefix and
        '000425' means that 'P000424' was the highest-numbered primer name
        sarting with P in the list of available primer names

        Parameters
        ----------
        prefix
          The prefix for the primers name

        available_primers_names
          List of already-allocated primer names
        """
        max_index = 0
        for name in available_primers_names:
            match = re.match(prefix + '(\d+)', name)
            if match is None:
                continue
            index = int(match.groups()[0])
            if (index > max_index):
                max_index = index
        return prefix + str(max_index + 1).zfill(n_digits)

    def find_part_name_in_record(self, record, index):
        """Find a part where the sequence appears in the provided record.

        This is used to provide primers infos ("where does this come from ?").

        Parameters
        ----------
        sequence
          An ATGC string representing a sequence (of a primer)

        record
          A single Biopython records where to find the sequence. The parts
          should be marked by features with a qualifier ``{'part': part_name}``

        Returns
        -------
        part_name, index
          The name of the part, and position in the part, where it was found.
          If the sequence is found in different parts, only the first find
          is returned.
        """
        for f in record.features:
            if 'part' not in f.qualifiers:
                continue
            part_name = ''.join(f.qualifiers['part'])
            strand = f.location.strand
            start, end = sorted([int(e) for e in [f.location.start,
                                                  f.location.end]])
            if start <= index <= end:
                distance = (index - start) if (strand >= 0) else (end - index)
                return (part_name, distance)
        else:
            raise ValueError(
                'The given index %d does not correspond to any part in "%s"'
                % (index, record.id)
            )

    def find_subsequence_in_records(self, sequence, records):
        """Find a part where the sequence appears in the provided records.

        This is used to provide primers infos ("where does this come from ?").
        This will look for either the sequence or its reverse-complement.

        Parameters
        ----------
        sequence
          An ATGC string representing a sequence (of a primer)

        records
          A list of Biopython records where to find the sequence. The parts
          should be marked by features with a qualifier ``{'part': part_name}``

        Returns
        -------
        part_name, index
          The name of the part, and position in the part, where it was found.
          If the sequence is found in different parts, only the first find
          is returned.
        """
        for r in records:
            ind = r.seq.find(sequence)
            if ind > 0:
                part_name, index = self.find_part_name_in_record(r, ind)
                break
            else:
                ind = r.seq.reverse_complement().find(sequence)
                if ind > 0:
                    part_name, index = self.find_part_name_in_record(r, ind)
                    break
        else:
            raise ValueError('sequence not found in records')
        return (part_name, index)

    def name_subsequence_in_records(self, sequence, records, prefix='P'):
        """Write a table of primers with columns 'name', 'sequence', etc.

        The columns after 'sequence' are one column per primer metadata, such
        as 'available', 'infos', etc.

        Parameters
        ----------
        primers
          A list of primers, or list of list, as returned by ``select_primers``

        csv_path
          The path to a csv file to write to. If None, no file is written, only
          the pandas dataframe of the table is returned.

        Returns
        -------
        dataframe
          The pandas dataframe of the table.
        """
        for r in records:
            ind = r.seq.find(sequence)
            if ind > 0:
                part_name, index = self.find_part_name_in_record(r, ind)
                break
            else:
                ind = r.seq.reverse_complement().find(sequence)
                if ind > 0:
                    part_name, index = self.find_part_name_in_record(r, ind)
                    break
        else:
            raise ValueError('sequence not found in records')
        return "%s%s_%04d" % (prefix, part_name, index)

    def plot_coverage(self, records, selected_primers, pdf_path,
                      close_figures=True, sort_matches_by=()):
        """Plot the predicted sequencing coverage for each construct."""
        # PLOT THE PREDICTED SEQUENCING COVERAGE FOR EACH CONSTRUCT

        with PdfPages(pdf_path) as pdf:
            rec_prim = list(zip(records, selected_primers))
            self.logger(message='Plotting coverages...')
            for record, primers in self.logger.iter_bar(record=rec_prim):
                matches_set = simulate_sequencing(
                    record, primers=primers, read_range=self.read_range,
                    linear=record.__dict__.get('linear', True)
                )
                if matches_set is None:
                    raise ValueError("No matches !")
                if len(sort_matches_by):
                    matches_set.sort_matches(sort_matches_by)
                ax = matches_set.plot(plot_reference=True, plot_coverage=False)

                ax.set_title(record.id, fontsize=14, weight='bold')
                for x in self.compute_indices_to_cover(record):
                    ax.axvline(x, c='#fceddb', lw=2, ls='--', zorder=-2000)
                pdf.savefig(ax.figure, bbox_inches='tight')
                if close_figures:
                    plt.close(ax.figure)

    def compute_coverage_records(self, records, selected_primers,
                                 target_dir=None):
        result = [
            simulate_sequencing(
                record, primers=primers, read_range=self.read_range,
                linear=record.__dict__.get('linear', True)
            ).to_biopython_record(record_id=record.id + '_coverage')
            for record, primers in zip(records, selected_primers)
        ]
        if target_dir is not None:
            root = flametree.file_tree(target_dir, replace=False)
            for record in result:
                SeqIO.write(record, root._file('%s.gb' % record.id), 'genbank')



    def write_records_primers_table(self, selected_primers, records, sep='|',
                                    csv_path=None):
        """Write a table with columns 'construct','primers' for this construct.

        Parameters
        ----------
        selected_primers
          The list of list of primers, as returned by the ``select_primers``
          method.

        records
          The list of records, as provided to the ``select_primers`` method.

        sep
          The separator between the different primer names in the ``primers``
          column. Avoid ';' or ',' as this might be used by the CSV formatter.

        csv_path
          The path to a csv file to write to. If None, no file is written, only
          the pandas dataframe of the table is returned.

        Returns
        -------
        dataframe
          The pandas dataframe of the table.
        """
        df = pandas.DataFrame.from_records([
            {
                'record': record.id,
                'primers': sep.join([p.name for p in primers_list])
            }
            for record, primers_list in zip(records, selected_primers)
        ], columns=['record', 'primers'])
        if csv_path is not None:
            df.to_csv(csv_path, index=False)
        return df

    def write_primers_table(self, selected_primers, csv_path=None):
        """Write a table of primers with columns 'name', 'sequence', etc.

        The columns after 'sequence' are one column per primer metadata, such
        as 'available', 'infos', etc.

        Parameters
        ----------
        primers
          A list of primers, or list of list, as returned by ``select_primers``

        csv_path
          The path to a csv file to write to. If None, no file is written, only
          the pandas dataframe of the table is returned.

        Returns
        -------
        dataframe
          The pandas dataframe of the table.
        """
        if isinstance(selected_primers[0], (list, tuple)):
            selected_primers = set([p for pp in selected_primers for p in pp])
        df = Primer.list_to_spreadsheet(selected_primers)
        df = df.sort_values(by=['available', 'name'])
        if csv_path is not None:
            df.to_csv(csv_path, index=False)
        return df

    def write_multifile_report(self, records, selected_primers,
                               sort_matches_by=('center', 'strand'),
                               target='@memory'):
        """Plot a full report in a folder or zip or in memory.

        The report contains:
        - A PDF 'coverages_plots.pdf', where each page shows how one construct
          should be covered by the primers' respective sequencing.
        - A spreadsheet 'primers_list.csv' indicating for each primer its
          sequence and whether it is already available.
        - A spreadsheet indicating, for each construct,

        Parameters
        -----------

        records
          A list of construct records

        selected_primers
          A list of list of primers (one list for each record)

        target
          Either the path to a directory or zip, or "@memory" to return a
          bytestring of binary data representing the zip file.

        """
        root = flametree.file_tree(target)
        self.plot_coverage(
            records=records,
            selected_primers=selected_primers,
            pdf_path=root._file('coverages_plots.pdf').open('wb'),
            sort_matches_by=sort_matches_by)
        self.write_primers_table(
            selected_primers=selected_primers,
            csv_path=root._file('primers_list.csv').open('w')
        )
        self.write_records_primers_table(
            records=records,
            selected_primers=selected_primers,
            csv_path=root._file('primers_per_record.csv').open('w')
        )
        return root._close()
