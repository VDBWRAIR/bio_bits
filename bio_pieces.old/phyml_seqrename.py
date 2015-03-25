import os
from os.path import *
import sys
import re
import argparse
from itertools import izip

'''
Garbage Test stuff
'''
test_infile = '1.H3N2_ConcatAll.fas'
test_phyml = 'input.phy'
test_tree = 'input.phy_phyml_tree.txt'

def main( ):
    args = parse_args()
    outstream = sys.stdout
    for renamefile in args.renamefile:
        rename_contents(
            args.fasta,
            renamefile,
            args.renamestr,
            outstream
        )

def numsortkey( item ):
    return int( re.search('\d+', item).group(0) )

def rename_contents( inputfile, renamefile, renamestr, outputstream=sys.stdout ):
    '''
        Rename the contents of a file
    '''
    rl = rename_list(inputfile, renamefile, renamestr)
    with open(renamefile) as fh:
        for line in fh:
            for find, replace in rl:
                line = line.replace(find, replace)
            outputstream.write( line )

def rename_list( renamefrom, renamefile, renamestr ):
    '''
    Reads renamefrom(fasta file) and pulls out all sequence identifiers.
    Then finds all renamestr patterns in renamefile

    Then returns a (find, replace) list from zipping those two lists

    Assumes that renamefrom and renamefile have the same number of
    items to zip together
    '''
    rl = get_rename_list( renamefile, renamestr )
    seqenum = get_seq_enumeration( renamefrom )

    if len(rl) == 0:
        raise Exception(
            'rename pattern \'{0}\' did not match anything in {1}'.format(
                renamestr, renamefile
            )
        )

    if len(rl) != len(seqenum):
        raise Exception(
            "There are more items to rename than there are available"
        )
    
    lst = []
    for rename, seq_num in izip(rl,seqenum):
        lst.append( (rename, seq_num[1]) )
    return lst

def get_rename_list( renamefile, renamestr ):
    '''
    Get a list of all renameable strings from a file. Should return a list
    of the same length as get_seq_enumeration from the fasta file
    Uses sortkey as the way to sort the list to ensure it is in order
    '''
    contents = ''
    with open(renamefile) as fh:
        contents = fh.read()
    rl = re.findall( renamestr, contents )
    return rl

def get_seq_enumeration( fastafile, start=0 ):
    '''
    Return enumeration of fasta sequence names
    '''
    seq_enumeration = []
    with open(fastafile) as fh:
        i = 0
        for line in fh:
            # Seq id starts with >
            if line.startswith('>'):
                # Only seqid without >
                seqid = line.strip()[1:]
                seq_enumeration.append( (i, seqid) )
                i += 1
    return seq_enumeration

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Rename Sequence names in another file ' \
            'based on an input fasta file'
    )

    parser.add_argument(
        '--renamestr',
        dest='renamestr',
        default='Seq\d+_',
        help='The string pattern to replace in the renamefile\'s[Default: %(default)s'
    )

    parser.add_argument(
        'fasta',
        help='Input fasta file in order to rename sequences'
    )

    parser.add_argument(
        'renamefile',
        nargs='+',
        help='List of files to have Seq##_ renamed in'
    )

    return parser.parse_args()
