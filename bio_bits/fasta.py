from __future__ import print_function
import sys
import argparse
import string

from Bio.SeqIO import parse
from Bio.SeqIO import FastaIO

def parse_args():
    parser = argparse.ArgumentParser(
        description='''Simplistic fasta manipulator'''
    )
    parser.add_argument(
        'fasta',
        help='Fasta file path or - to read from standard input'
    )
    parser.add_argument(
        '--wrap', default=False, action='store_const', const=80,
        help='Default action. Converts all multi line sequences to single lines'
    )
    parser.add_argument(
        '--split', action='store_true', default=False,
        help='If set then split the input fasta on each identifier '
            ' and create a new file for each identifier'
    )
    return parser.parse_args()

def writefasta(infile, outfile=sys.stdout, wrap=None):
    fasta_out = FastaIO.FastaWriter(outfile, wrap=wrap)
    fasta_out.write_file(parse(infile, 'fasta'))

def splitfasta(infile, wrap=False):
    for seq in parse(infile, 'fasta'):
        outfile = str(seq.id)
        for p in string.punctuation:
            outfile = outfile.replace(p, '_')
        outfile = outfile + '.fasta'
        with open(outfile, 'w') as fh:
            fasta_out = FastaIO.FastaWriter(fh, wrap=wrap)
            fasta_out.write_header() # Does nothing, but required
            fasta_out.write_record(seq)
            fasta_out.write_footer() # Does nothing, but required

def main():
    args = parse_args()
    in_fasta = args.fasta

    # Normalize input stream
    if in_fasta == '-':
        _input = sys.stdin
    else:
        _input = open(in_fasta)

    if args.split:
        splitfasta(_input, args.wrap)
    else:
        writefasta(_input, sys.stdout, wrap=args.wrap)
