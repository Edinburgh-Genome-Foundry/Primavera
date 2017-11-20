from primavera import PrimerSelector, Primer, load_record
import os

# LOAD ALL RECORDS TO ANALYSE AND AVAILABLE PRIMERS
constructs_path = os.path.join('data', 'constructs')
records = [
    load_record(os.path.join(constructs_path, f), linear=False)
    for f in sorted(os.listdir(constructs_path)) if f.endswith('.gb')
]
primers_path = os.path.join('data', "available_primers.fa")
available_primers = Primer.list_from_fasta(primers_path)


# SELECT THE BEST PRIMERS
selector = PrimerSelector(read_range=(150, 800), tm_range=(55, 70))
selected_primers = selector.select_primers(records, available_primers)


# PLOT THE PREDICTED SEQUENCING COVERAGE FOR EACH CONSTRUCT
selector.plot_coverage(records=records, selected_primers=selected_primers,
                       pdf_path='primer_selection_example.pdf')


# WRITE ALL PRIMERS IN A CSV FILE (PRIMERS TO ORDER ARE FIRST)
selector.write_primers_table(selected_primers=selected_primers,
                             csv_path='primer_selection_example.csv')
