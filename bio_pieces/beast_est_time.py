#!/bin/bash

# Estimates the amount of time remaining in a beast run
# You have to get the number of chains from the xml file
# Let's say that your numchains is 1,000,000
# You would run as follows
#
# $> python beast_est_time.py 1000000 <paste last beast line here>

import re
import sys

def hours_states_to_sec(line, chainlength):
    r'''
    >>> line = '20000   -61227.3011     -30325.1333     -30902.1678     932.046         1.60076E-4      0.35573         3487.00         3.39 hours/million states'
    >>> chainlength = 500000000
    >>> hours_states_to_sec(line, chainlength)
    6101755
    >>> hours_states_to_sec(line.replace('million','billion'), chainlength*1000)
    6101999
    >>> line = '20000   -61227.3011     -30325.1333     -30902.1678     932.046         1.60076E-4      0.35573         3487.00         0.39 hours/million states\n'
    >>> hours_states_to_sec(line, chainlength)
    '''
    p = '^(\d+).*?((\d+\.\d+)\shours/[mb]illion\sstates)$'
    m = re.search(p, line)
    if not m:
        raise ValueError(line + ' is not a valid beast output line')

    completed_states, stateline, states = m.groups()
    states = float(states)
    completed_states = float(completed_states)

    multiplier = 1000000
    if 'billion' in stateline:
        multiplier = multiplier * 1000

    states_to_complete = (float(chainlength) - completed_states) / multiplier
    sec = states * states_to_complete * 3600
    return int(sec) 

def sec_to_time(secs):
    '''
    >>> sec_to_time(1800)
    '0d 00:30:00'
    >>> sec_to_time(3600)
    '0d 01:00:00'
    >>> sec_to_time(86400)
    '1d 00:00:00'
    >>> sec_to_time(86400+3600+1800+60+10)
    '1d 01:31:10'
    '''
    days, hours = divmod(secs, 86400)
    hours, minutes = divmod(hours, 3600)
    minutes, seconds = divmod(minutes, 60)
    return "{0}d {1:02}:{2:02}:{3:02}".format(days,hours,minutes,seconds)

def main():
    if len(sys.argv) > 1:
        states = sys.argv[1]
        line = ' '.join(sys.argv[2:])
        s = hours_states_to_sec(line, states)
        print sec_to_time(s)
    else:
        print "Usage: beast_est_time.py <numchains> <beastline>"

if __name__ == '__main__':
    main()
