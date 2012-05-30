#!/usr/bin/env python2.4

##########################################################################
##                       Sequence Concat
##	Author: Tyghe Vallard						                        
##	Date: 5/30/2012							                            
##	Version: 1.0							                            
##	Description:							                            
##      This script will merge all sequences with identical names in a
##      file.
##
##  Example:
##      SampleFile1.fasta contents
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##      >identifier1
##      TTTTTTTTTTTTTTTTTTTTTTTTT
##      >identifier1
##      GGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCC
##      >identifier2
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##
##      OUTPUT
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAA
##      
##########################################################################
import os
import optparse import OptionParser
import sys

#==================== Variables =========================
#---------------------------------------------------------

#==================== Functions =========================
def set_opts( ):
    """ Set script options """
    global parser
    global options

    parser.add_option( "-d", "--fasta_dir", dest="fasta_dir", help="The fasta directory of the sequences to merge" )
    parser.add_option( "--debug", action="store", dest="debug", default=0, help="Set debug level. Normal(default) = 0, Warn = 1, Info = 2, Debug = 3" )
    parser.add_option( "--strip", dest="strip_chars", default="-", help="List of characters to strip from the sequences. Default is '-'" )
    options,args = parser.parse_args()
    if not options.fasta_dir:
        parser.print_help()
        parser.error( "Need to specify the fasta directory" )

def main( argv ):
    set_opts()
#---------------------------------------------------------

if __name__ == '__main__':
    main( argv )
