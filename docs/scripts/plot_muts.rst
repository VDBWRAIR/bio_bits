plot_muts
=========

Usage
-----

You can view the usage of degen_regions via::

    plot_muts --help
    
Example
-------

.. code-block:: bash

    plot_muts --refs tests/testinput/refs.fas --query tests/testinput/query.fas --out plot.png

The ``--out`` option is optional. If it is not provided, the plot will pop up on 
the user's screen automatically. If this does not work, try saving the image using ``--out`` instead.

Input File Requirements
-----------------------

The input must be fasta format. Both the query and ref files can have any number of sequences.

The year should be the last part of the ID, preceded by an underscore. e.g.::

    >some|info|blah_blah_1995
    >some_1995
    
If the ID uses '/' rather than underscore, plot_muts currently accepts the year 
as the *fourth* field. e.g.::

    >some/info/blah/1995
    >some/info/blah/1995/more/info
