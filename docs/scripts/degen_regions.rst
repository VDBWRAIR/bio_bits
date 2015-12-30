degen_regions
=============

Finds all degenerate bases in a given fasta input file that may contain multiple
sequeces and reports their position as well as the annotated gene name that contains
them.

The fasta file must be previously aligned to the query sequence. That is, if you
are using a genbank annotation file or having the script download it for you, you
should have aligned all your input sequences to that sequence.

The annotation is retrieved via supplied genbank accession, genbank file path or
gene tab/csv file.

Usage
-----

You can view the usage of degen_regions via::

    degen_regions --help

Using Genbank Files
-------------------

If you already have downloaded the genbank annotation file(typically the extension
is .gb) you can use the `--gb-file` argument

The following will use the test input fasta file as well as the test input
genbank file to find all degenerate bases and will put the output in a tab separated
file called output.tsv

.. code-block:: bash

    degen_regions -i tests/Den4_MAAPS_TestData16.fasta -o output.tsv --gb-file tests/testinput/sequence.gb

Fetching Genbank Files Automatically
------------------------------------

If you want the script to automatically fetch the Genbank annotation file from the
internet you can use the `--gb-id` option and specify an accession number.

.. code-block:: bash

    degen_regions -i tests/Den4_MAAPS_TestData16.fasta -o output.tsv --gb-id KJ189367

Using tab/csv file of gene annotation info
------------------------------------------

If you have a tab/csv file of gene annotations you can supply that using the
`--tab-file` argument

You can read more about the format of the tab/csv annotation file in the :doc:`degen` docs

.. code-block:: bash

    degen_regions -i tests/Den4_MAAPS_TestData16.fasta -o output.tsv --gb-file tests/testinput/sequence.gb

Manually specify CDS
--------------------

You can use the ``--cds`` argument to set the coding region.
This argument should be comma separated such as ``start,stop``.
Specifying this argument will override any other cds found in the tab file, genbank
file or fetched genbank file.

The following would mark all locations as NON-CODING as you are specifying that only
position 1 is coding

.. code-block:: bash

    degen_regions -i tests/Den4_MAAPS_TestData16.fasta -o output.tsv --gb-file tests/testinput/sequence.gb --cds 1,1

Output
------

The output is a simple tab separated file

.. include:: ../../tests/testinput/ctl_expected.tsv
    :literal:
    :end-before: 1909_Den
