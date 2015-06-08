beast_wrapper
=============

Beast wrapper is intended as a helper script to run beast. At this point it just
runs beast with the same arguments you would normally give to beast from the command
line and just adds a estimated time left column to the console output

Example
-------

.. code-block:: bash

    $> beast_wrapper -beagle_SSE my_beast.xml
    ...
    state   Posterior       Prior           Likelihood      rootHeight      Den1_RO1_genomes_03Feb15_aln.ucld.mean  location.clock.rate location.nonZeroRates
    0   -86527.5880     -6850.8316      -79676.7564     57.6772         1.16103E-3      4.86012         15.0000         -
    20000   -29044.3753     -1123.5287      -27920.8466     288.102         3.02471E-4      0.11891         16.0000         0.21 hours/million states   2d 04:29:44
    40000   -25517.9525     -979.5343       -24538.4182     211.705         1.35118E-4      0.25060         16.0000         0.25 hours/million states   2d 14:29:24
    60000   -24212.1250     -1040.4103      -23171.7147     188.454         1.05572E-4      0.18908         15.0000         0.25 hours/million states   2d 14:29:06
    80000   -24097.9354     -1019.8099      -23078.1256     182.242         1.53593E-4      0.12857         16.0000         0.26 hours/million states   2d 16:58:45
    100000  -24121.5382     -1105.6545      -23015.8837     178.060         1.26907E-4      0.10367         17.0000         0.27 hours/million states   2d 19:28:22
    120000  -23930.6897     -1105.7390      -22824.9507     187.411         1.01885E-4      0.34214         17.0000         0.27 hours/million states   2d 19:28:03
    140000  -23869.4856     -1087.1915      -22782.2942     178.535         8.76375E-5      0.26128         18.0000         0.26 hours/million states   2d 16:57:48
