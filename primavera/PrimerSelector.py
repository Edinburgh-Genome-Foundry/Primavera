from copy import deepcopy
import re
from collections import defaultdict

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

from proglog import TqdmProgressBarLogger, ProgressBarLogger
from dnachisel import (AvoidPattern, repeated_kmers, homopolymer_pattern,
                       DnaOptimizationProblem)

from .sequencing_simulation import simulate_sequencing
from .biotools import (reverse_complement, find_non_unique_segments,
                       find_best_primer_locations)
from .tools import minimal_cover, segments_to_array, group_overlapping_segments
from .Primer import Primer

class PrimerSelectorLogger(TqdmProgressBarLogger):

    def __init__(self, bars=('record', 'primer'), notebook='default'):
        ignored_bars = set(('record', 'primer')).difference(bars)
        TqdmProgressBarLogger.__init__(self, bars=bars, notebook=notebook,
                                       ignored_bars=ignored_bars,
                                       min_time_interval=0.2)

class PrimerSelector:

    def __init__(self, read_range=(150, 800), size_range=(16, 25),
                 tm_range=(55, 70), primer_conditions=(),
                 primer_reuse_bonus=2, logger='bars',
                 coverage_resolution=5):
        self.read_range = read_range
        self.size_range = size_range
        self.tm_range = tm_range
        self.primers_conditions = primer_conditions
        self.coverage_resolution = coverage_resolution
        self.primer_reuse_bonus = 2
        if logger == 'bars':
            logger = PrimerSelectorLogger()
        if logger is None:
            logger = ProgressBarLogger()
        self.logger = logger

    def select_primers(self, records, available_primers=(),
                       new_primers_prefix='P', new_primers_digits=6):

        available_primers_dict = {p.sequence: p for p in available_primers}
        available_primers_seqs = set([p.sequence for p in available_primers])

        # COMPUTE PRIMERS AND COVERAGES
        indices_to_cover = {}
        primers_coverages = defaultdict(lambda *a: set())
        self.logger(message='Analyzing the records...')
        for record in self.logger.iter_bar(record=records):
            indices_to_cover[record.id] = {
                ind: '%s_%03d' % (record.id, ind)
                for ind in self.compute_indices_to_cover(record)
            }
            coverages = self.compute_all_primers_coverage_on_record(
                record, available_primers=available_primers_seqs,
                indices_to_cover=indices_to_cover[record.id])
            for primer, coverage in coverages.items():
                primers_coverages[primer].update(coverage)

        # FIND GLOBAL MINIMAL COVER
        self.logger(message='Selecting primers...')
        elements_set = set(
            index
            for rec_id, named_indices in indices_to_cover.items()
            for index in named_indices.values()
        )
        def heuristic(named_subset, selected):
            name, subset = named_subset
            primer_is_reused = name in available_primers_seqs
            reuse_bonus = self.primer_reuse_bonus * primer_is_reused
            return len(subset) + reuse_bonus
        #subsets = deepcopy(list(primers_coverages.items()))
        subsets = primers_coverages.items()
        primers_cover = minimal_cover(elements_set, subsets=subsets,
                                      heuristic=heuristic)

        # REORGANIZE AND NAME THE SELECTED PRIMERS
        available_primers_names = [p.name for p in available_primers]
        selected_primers = []
        selected_primer_from_seq = {}
        for primer_seq in primers_cover:
            if primer_seq in available_primers_dict:
                name = available_primers_dict[primer_seq].name
                primer = available_primers_dict[primer_seq]
                primer = Primer(primer.name, primer.sequence,
                                metadata={'available': True})
            else:
                name = self.generate_primer_name(
                    available_primers_names=available_primers_names,
                    prefix=new_primers_prefix,
                    n_digits=new_primers_digits
                )
                available_primers_names.append(name)
                primer = Primer(name, primer_seq,
                                metadata={'available': False})
            selected_primers.append(primer)
            selected_primer_from_seq[primer_seq] = primer

        # CHOOSE A MINIMAL PRIMER COVER FOR EACH CONSTRUCT
        per_record_primers = []
        for record in self.logger.iter_bar(record=records):
            elements = set(indices_to_cover[record.id].values())
            subcovers = {
                prim_seq: primers_coverages[prim_seq].intersection(elements)
                for prim_seq in primers_cover
            }
            # subsets = deepcopy(list(subcovers.items()))
            subsets = list(subcovers.items())
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
        segments_to_cover = [
            sorted([int(f.location.start), int(f.location.end)])
            for f in record.features
            if f.type == 'misc_feature'
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
        pattern_constraints = [AvoidPattern(homopolymer_pattern(c, 5))
                               for c in 'ATGC']
        kmer_constraints = [AvoidPattern(repeated_kmers(k, n))
                            for k, n in [(4, 2), (3, 3), (2, 4)]]
        problem = DnaOptimizationProblem(
            sequence=record,
            constraints=pattern_constraints + kmer_constraints
        )
        constraints_breaches = group_overlapping_segments([
            (f.location.start, f.location.end)
            for f in problem.constraints_evaluations().locations_as_features()
        ])
        return segments_to_array(constraints_breaches, len(record))

    def compute_user_forbidden_locations(self, record):
        forbidden_segments = [
            sorted([int(f.location.start), int(f.location.end)])
            for f in record.features
            if f.type == 'misc_feature'
            and "".join(f.qualifiers.get('label', '')) == 'no_primer'
        ]
        return segments_to_array(forbidden_segments, len(record))

    def compute_nonunique_segments_locations(self, record):
        sequence = str(record.seq)
        non_unique_segments = find_non_unique_segments(sequence)
        return segments_to_array(non_unique_segments, len(record))

    def compute_all_forbidden_locations(self, record):
        return np.maximum(*(f(record) for f in (
            self.compute_nonunique_segments_locations,
            self.compute_forbidden_patterns_locations,
            self.compute_user_forbidden_locations
        )))

    def compute_sequence_primers(self, record):
        sequence = str(record.seq)
        L = len(sequence)
        rev_sequence = reverse_complement(sequence)
        locations = find_best_primer_locations(
            sequence,
            tm_range=self.tm_range,
            size_range=self.size_range
        )
        return list(set(sum([
            [
                (sequence[l[0]: l[1]], (l[0], l[1], 1)),
                (rev_sequence[L - l[1]: L - l[0]], (l[0], l[1], -1))
            ]
            for l in locations
            if l is not None
        ], [])))

    def compute_all_primers_coverage_on_record(self, record, indices_to_cover,
                                               available_primers):
        sequence = str(record.seq)
        linear = record.__dict__.get('linear', True)
        L = len(sequence)
        forbidden_locations = self.compute_all_forbidden_locations(record)
        sequence_primers = self.compute_sequence_primers(record)
        reusable_primers_sequences = [
            (primer, location) for (primer, location) in [
                (primer, self.locate_primer_sequence(primer, sequence))
                for primer in available_primers
            ]
            if location is not None
        ]
        primers = reusable_primers_sequences + sequence_primers
        iter_bar = self.logger.iter_bar
        primers_coverages = {}
        for primer, (start, end, strand) in iter_bar(primer=primers):
            if forbidden_locations[start:end].max():
                continue
            if not all(cond(primer) for cond in self.primers_conditions):
                continue

            if strand == 1:
                coverage_start = end + self.read_range[0]
                coverage_end = end + self.read_range[1]
            else:
                coverage_start = start - self.read_range[1]
                coverage_end = start - self.read_range[0]
            if (not linear) and (coverage_start < 0):
                a, b, c, d = -1000, coverage_end, L + coverage_start, L + 1000
            elif (not linear) and (coverage_end > L):
                a, b, c, d = -1000, coverage_end - L, coverage_start, L + 1000
            else:
                a, b, c, d = (coverage_start, coverage_end,
                              coverage_start, coverage_end)
            primers_coverages[primer] = set([
                indices_to_cover[ind]
                for ind in indices_to_cover
                if (a <= ind <= b) or (c <= ind <= d)
            ])
        return primers_coverages

    def find_part_name_in_record(self, record, index):
        for f in record.features:
            if 'part' not in str(f.qualifiers):
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

    def generate_primer_name(self, prefix='P', available_primers_names=(),
                             n_digits=6):
        max_index = 0
        for name in available_primers_names:
            if not name.startswith(prefix):
                continue
            index = int(re.search(r'\d+', name[len(prefix):]).group())
            if (index is not None) and (index > max_index):
                max_index = index
        return prefix + str(max_index + 1).zfill(n_digits)


    def name_subsequence_in_records(self, sequence, records, prefix='P'):
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

    def plot_coverage(self, records, selected_primers, pdf_path):
        """Plot the predicted sequencing coverage for each construct."""
        # PLOT THE PREDICTED SEQUENCING COVERAGE FOR EACH CONSTRUCT

        with PdfPages(pdf_path) as pdf:
            iterator = zip(records, selected_primers)
            self.logger(message='Plotting coverages...')
            for record, primers in self.logger.iter_bar(record=iterator):
                matches_set = simulate_sequencing(
                    record, primers=primers, read_range=self.read_range,
                    linear=record.__dict__.get('linear', True)
                )
                ax = matches_set.plot(plot_reference=True, plot_coverage=False)
                ax.set_title(record.id, fontsize=14, weight='bold')
                for x in self.compute_indices_to_cover(record):
                    ax.axvline(x, c='#fceddb', lw=2, ls='--', zorder=-2000)
                pdf.savefig(ax.figure, bbox_inches='tight')
                plt.close(ax.figure)

    def write_primers_table(self, primers, csv_path=None):
        if isinstance(primers[0], (list, tuple)):
            primers = set([p for pp in primers for p in pp])
        df = Primer.list_to_spreadsheet(primers)
        df = df.sort_values(by=['available', 'name'])
        if csv_path is not None:
            df.to_csv(csv_path, index=False)
        return df
