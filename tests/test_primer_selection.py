import os
import matplotlib
matplotlib.use("Agg")
from primavera import (load_record, Primer, PrimerSelector)


constructs_path = os.path.join('tests', 'data', 'constructs')
primers_path = os.path.join('tests', 'data', "available_primers.fa")

def test_PrimerSelector(tmpdir):
    records = [
        load_record(os.path.join(constructs_path, f), linear=False)
        for f in sorted(os.listdir(constructs_path)) if f.endswith('.gb')
    ][:3]
    available_primers = Primer.list_from_fasta(primers_path)

    # SELECT THE BEST PRIMERS
    selector = PrimerSelector(read_range=(150, 800), tm_range=(55, 70))
    selected_primers = selector.select_primers(records, available_primers)

    # WRITE ALL PRIMERS IN A CSV FILE (PRIMERS TO ORDER ARE FIRST)
    csv_path = os.path.join(str(tmpdir), 'primer_selection_example.csv')
    table = selector.write_primers_table(selected_primers=selected_primers,
                                         csv_path=csv_path)
    assert (table.available.sum() == 2)

    selector.write_records_primers_table(selected_primers=selected_primers,
                                         records=records, csv_path=csv_path)

    # PLOT THE PREDICTED SEQUENCING COVERAGE FOR EACH CONSTRUCT
    pdf_path = os.path.join(str(tmpdir), 'primer_selection_example.pdf')
    selector.plot_coverage(records=records, selected_primers=selected_primers,
                           pdf_path=pdf_path)
    selector.write_multifile_report(records, selected_primers,
                                    target=str(tmpdir))
