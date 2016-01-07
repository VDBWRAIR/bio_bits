fasta
=====

fasta is a very simple script to help mangle fasta files. Currently it only supports
the ability to convert fasta files that have sequences that span multiple new lines
into single lines.

Later on, it may be expanded more to include even more useful fasta features.

Usage
-----

.. code-block:: bash

    fasta --help

Examples
^^^^^^^^

The following examples all use the test fasta file found under 
``tests/testinput/col.fasta``

.. include:: ../../tests/testinput/col.fasta
    :literal:

Read fasta input from standard input
++++++++++++++++++++++++++++++++++++

The following could output the fasta sequences as one line to your terminal(stdout)
but reading from the pipe. This is useful if you want to use it in a pipeline.

.. code-block:: bash

    $> cat tests/testinput/col.fasta | fasta -

Read fasta input from file
++++++++++++++++++++++++++

The following could output the fasta sequences as one line to your terminal(stdout)
as well, but reading straight from the file.

.. code-block:: bash

    $> fasta tests/testinput/col.fasta

Simple shell pipeline using fasta
+++++++++++++++++++++++++++++++++

The following is a simple shell pipeline to count how many A's there are in the
sequence lines. There should be 160 since ``col.fasta`` is 80 characters per line
and only the first line of each sequence has ``A`` and there are 2 sequences.

.. code-block:: bash

    $> fasta tests/testinput/col.fasta | grep -v '>' | grep -Eo '[Aa]' | wc -l
    160
