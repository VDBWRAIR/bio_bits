fasta
=====

fasta is a very simple script to help mangle fasta files. 

* Supports converting multiline sequences into single line
* Supports splitting fasta file into separate files each named after the identifier
* Supports disambiguating ambiguous sequences

Usage
-----

.. code-block:: bash

    fasta --help

Examples
--------

The following examples all use the test fasta file found under 
``tests/testinput/col.fasta``

.. include:: ../../tests/testinput/col.fasta
    :literal:

Convert column fasta into single lines
++++++++++++++++++++++++++++++++++++++

The following could output the fasta sequences as one line to your terminal(stdout)
but reading from the pipe. This is useful if you want to use it in a pipeline.

.. code-block:: bash

    $> cat tests/testinput/col.fasta | fasta -

The following could output the fasta sequences as one line to your terminal(stdout)
as well, but reading straight from the file.

.. code-block:: bash

    $> fasta tests/testinput/col.fasta

Convert single line fasta into column fasta
+++++++++++++++++++++++++++++++++++++++++++

The following would convert single line fasta sequences into column formatted
fasta. It defaults to using 80 characters for each column

.. code-block:: bash

    $> fasta tests/testinput/col.fasta

You can verify that it is wrapping correctly by simply piping the fasta command
back into itself from ``convert-to-single | wrap`` and then comparing to the original
input file.

.. code-block:: bash

    $> cat tests/testinput/col.fasta | fasta - | fasta --wrap - > newfile.fasta
    $> diff tests/testinput/col.fasta newline.fasta

There will be no output as there is no difference between ``newfile.fasta`` and
``tests/testinput/col.fasta``

Simple shell pipeline using fasta
+++++++++++++++++++++++++++++++++

The following is a simple shell pipeline to count how many A's there are in the
sequence lines. There should be 160 since ``col.fasta`` is 80 characters per line
and only the first line of each sequence has ``A`` and there are 2 sequences.

.. code-block:: bash

    $> fasta tests/testinput/col.fasta | grep -v '>' | grep -Eo '[Aa]' | wc -l
    160

Split fasta file into separate files named after identifiers
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following example shows how you can split a fasta file into multiple fasta
files each named after an identifier in the original

.. code-block:: bash

    $> fasta tests/testinput/col.fasta --split
    $> ls -l *.fasta
    sequence1.fasta
    sequence2____________________________.fasta

*Note* The reason sequence2 has such a long name is because it is replacing
all punctionation characters with underscores. ``col.fasta`` is a test file
that has a bunch of punctuation, hence all the underscores.

Similar to above, you can use input from standard input as the fasta input file

.. code-block:: bash

    $> cat tests/testinput/col.fasta | fasta --split -
    $> ls -l *.fasta
    sequence1.fasta
    sequence2____________________________.fasta

Disambiguate ambiguous sequences
++++++++++++++++++++++++++++++++

You can turn sequences that have ambiguous bases in them into all permutations
of the same sequence with the ambiguous bases turned into non-ambiguous bases.

There is an upper limit of 100 for how many sequences can be generated to avoid
creating thousands of sequences or consuming all of your computer's RAM.

If a sequence would generate more than 100 sequences, it will generate a message
such as::

    Sequence too_many has 7 ambiguous bases that would produce 128 permutations and was skipped

and it will be skipped.

.. code-block:: bash

    $> fasta --disambiguate tests/testinput/ambiguous.fasta > disambiguous.fasta
