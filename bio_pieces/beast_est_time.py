#!/bin/bash

# Estimates the amount of time remaining in a beast run
# You have to get the number of chains from the xml file
# Let's say that your numchains is 1,000,000
# You would run as follows
#
# $> python beast_est_time.py 1000000 <paste last beast line here>

import re
import sys


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
