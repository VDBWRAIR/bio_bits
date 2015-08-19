group_references
================

group_references splits an alignment file by reference into seperate FASTQ files. group_references takes a SAM or BAM file as input, and can optionally be given an output directory where the FASTQ files will be saved. If not output directory name is provided, the files will be saved in the new folder group_references_out.

.. code-block:: bash

    $> group_references contigs.bam 
    $> group_references contigs.bam --outdir split_fastqs
