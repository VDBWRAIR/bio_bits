beast_checkpoint
================

beast_checkpoint is a fork of https://gist.github.com/trvrb/5277297 that has been
rewritten in python and slightly improved as the ruby script seemed to have a few
errors.

It accepts any previously run or terminated beast run and will generate an xml
file that essentially starts from the last generated tree/log state.

Since beast is random in nature, there does not appear to be a way to restart the run
exactly from the same state that it left off.

Example
-------

We will use the benchmark2.xml file that comes with Beast 1.8
This file is located in::

    BEASTv1.8.0/examples/Benchmarks/benchmark2.xml


First you need to fix the benchmark2.xml because each taxa has a trailing space and
that is annoying

.. code-block:: bash

    $> sed 's/ "/"/' benchmark2.xml > beast.xml

Now run beast for about half of the iterations and hit CTRL-C to kill it
This benchmark is set to run 1,000,000 iterations so around 500,000 you can kill it.
Notice we are using a predefined seed

.. code-block:: bash

    $> seed=1234567890
    $> mkdir run1
    $> cp beast.xml run1/beast.xml
    $> beast -seed $seed -beagle_SSE beast.xml

Now we will want to re-run beast from that last state. We can use beast_checkpoint
to do so by supplying the original xml and the produced trees and log files.
We will put the new xml into a new directory since the .trees and .log files would
create an error or possibly be overwritten.

*NOTE* 
If your fileLog and treeFileLog do not have the same logEvery then when beast
exits you may end up with more/less tree states than log states. 
For now you will have to manually edit the files and ensure that the last tree state
matches the last log state.

.. todo::

    Could be possible to get beast_checkpoint to check for that scenario and 
    use the last tree state that matches the last log state

.. code-block:: bash

    $> mkdir run2
    $> beast_checkpoint beast.xml *.trees *.log > run2/beast.xml

Now you can simply just re-run beast on the new xml using the same seed

.. code-block:: bash

    $> cd run2
    $> beast_checkpoint -seed $seed -beagle_SSE beast.xml

Tracer
------

If you name your runs sequentially as we did in the example(aka, run1, run2,...)
then you can easily load all log files into tracer via the command line as follows

.. code-block:: bash

    tracer run*/*.log

LogCombiner
-----------

After you have run all your beast checkpointed xml files you will probably want to
combine them with logcombiner which comes with beast
