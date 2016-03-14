fasta
=====

fasta is a very simple script to help mangle fasta files. 

* Supports converting multiline sequences into single line
* Supports splitting fasta file into separate files each named after the identifier

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

The following is a simple shell pipeline using fasta to ensure all sequences are
on one line

.. code-block:: bash

    $> cat tests/testinput/col.fasta | fasta -

Or if you want to you can read straight from a fasta file

.. code-block:: bash

    $> fasta tests/testinput/col.fasta

Convert single line fasta into column fasta
+++++++++++++++++++++++++++++++++++++++++++

The following would convert single line fasta sequences into column formatted
fasta. It defaults to using 80 characters for each column

.. code-block:: bash

    $> fasta tests/testinput/col.fasta

You can verify that it is wrapping correctly by simply piping the fasta command
back into itself and then comparing to the original input file.

Here you can see we do that and then use diff to show there is no difference between
the original file(col.fasta) and the new one(newline.fasta)

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
