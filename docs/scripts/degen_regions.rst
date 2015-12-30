degen_regions
=============

Finds all degenerate bases in a given fasta input file that may contain multiple
sequeces and reports their position as well as the annotated gene name that contains
them.

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
`--tab_file` argument

You can read more about the format of the tab/csv annotation file in the :doc:`degen` docs

.. code-block:: bash

    degen_regions -i tests/Den4_MAAPS_TestData16.fasta -o output.tsv --gb-file tests/testinput/sequence.gb

Output
------

The output is a simple tab separated file

.. include:: ../../tests/testinput/ctl_expected.tsv
    :literal:
    :end-before: 948_Den4/AY618992_1/Thailand/2001/Den4_1
