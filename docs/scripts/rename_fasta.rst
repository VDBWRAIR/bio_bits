rename_fasta
============

Many times you find you have a fasta file where the identifiers are all wrong and you
want to rename them all via some mapping file.

Take the example where you have the following fasta file(example.fasta)::

    >id1
    ATGC
    >id2
    ATGC
    >id3
    ATGC

You want to rename each identifier(id1, id2, id3) based on a mapping you have.
In a file called renamelist.csv you would have the following::

    #From,To
    id1,samplename1
    id2,samplename2
    id3,samplename3

Then to rename your fasta without replacing the original file you have two options:

#. Rename without replacing original file

    .. code-block:: bash

        rename_fasta renamelist.csv example.fasta > renamedfasta.fasta

#. Rename replacing original file's contents

    .. code-block:: bash

        reanme_fasta renamelist.csv example.fasta --inplace

Rename Mapping File Syntax
--------------------------

The file you specify as the rename map file is a simple comma separated text file.

The following rules apply to the format:

* The first entry is the identifier to find in the supplied fasta file.
* The second entry is what to replace the found identifier with
* Any line beginning with a pound sign(#) will be ignored by the renamer

Missing identifiers that are in fasta but not rename file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case where your fasta file contains an identifier that is not in the rename
map file you supply, an error will be displayed in the console telling you as such::

    idwhatever is not in provided mapping
