.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Primavera/master/docs/_static/images/title.png
   :alt: [logo]
   :align: center
   :width: 600px

Primavera
==========

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/Primavera.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/Primavera
   :alt: Travis CI build status

Primavera is a Python library to plan and analyze primer-based verification of DNA assemblies, using Sanger sequencing or verification PCR. It implements methods to design and select primers to ensure that the relevant assembly segments will be covered, and it

Usage
-----

**For primer selection**, assuming

.. code:: python

    from primavera import PrimerSelector, Primer, load_record
    import os

    # LOAD ASSEMBLIES RECORDS AND AVAILABLE PRIMERS
    records = [load_record(file_path, linear=False)
               for file_path in ['my_record_1.gb', 'my_record_2.gb'...]]
    available_primers = Primer.list_from_fasta("example_primers.fa")

    # SELECT THE BEST PRIMERS
    selector = PrimerSelector(tm_range=(55, 70), size_range=(16, 25))
    selected_primers = selector.select_primers(records, available_primers)

    # PLOT THE COVERAGE AND WRITE THE PRIMERS IN A SPREADSHEET
    selector.plot_coverage(records, selected_primers, 'coverage.pdf')
    selector.write_primers_table(selected_primers, 'selected_primers.csv')


Installation
-------------

(Soon) You can install DnaChisel through PIP

.. code::

    sudo pip install primavera

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

License = MIT
--------------

Primavera is an open-source software originally written at the `Edinburgh Genome Foundry <http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_ and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

Contribute
-----------

Primavera is an open-source software originally developed at the Edinburgh
Genome Foundry by Zulko and released on Github under the MIT licence (copyright Edinburgh Genome Foundry). Everyone is welcome to contribute !
