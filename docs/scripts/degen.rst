degen
=====

Find genes where a sequence has degenerate bases.

How-to
------

Usage:
    degen.py <fasta> <options>

Options:

    --gb-id=<accession_id>   Accession id for reference
    --gb-file=<gbfile>       Local Genbank file for reference
    --tab-file=<tabfile>     TSV/CSV file for reference with fields name,start,end

Example:
--------

    .. code-block:: bash

        degen sequence.fasta --gb-id   12398.91
        degen sequence.fasta --gb-file  tests/testinput/sequence.gb
        degen sequence.fasta --tab-file tests/testinput/degen.tab
        degen sequence.fasta --tab-file tests/testinput/degen.csv


Output:
-------

Gene name, degnerate position, degenerate base::

    anchored capsid protein         85      R
    anchored capsid protein         88      Y
    membrane glycoprotein precursor 509     R
    nonstructural protein NS5       8513    Y
    nonstructural protein NS5       8514    Y
    nonstructural protein NS5       8515    Y
    anchored capsid protein         85      R
    anchored capsid protein         88      Y
    membrane glycoprotein precursor 509     R
    nonstructural protein NS5       8513    Y
    nonstructural protein NS5       8514    Y
    nonstructural protein NS5       8515    Y

Gene/Tab File
-------------

degen.tab could look like::

    genename    start    stop
    foo    1    2
    bar    9    33 

The headers do not matter, but the start field must always come before the stop 
field, so the below example would also be valid::

    start    GENEName    stop
    1    foo    2
    9    bar    33 

or optionally without headers::

    1    foo    2
    9    bar    33 

alternatively, with commas in place of tabs::

    name,start,stop
    foo,1,2
    bar,9,33 

You can also specify a coding region(CDS) in your file as well::

    name,start,stop

    CDS,3,33
    foo,1,2
    bar,9,33 

Genbank File
------------

As downloaded from NCBI's entrez database.
Use this option if you don't have internet access.

An example

.. include:: ../../tests/testinput/sequence.gb
    :literal:
