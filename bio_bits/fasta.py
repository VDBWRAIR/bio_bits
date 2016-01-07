from __future__ import print_function
import sys
import argparse

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
    return parser.parse_args()

def main():
    args = parse_args()
    in_fasta = args.fasta

    # Normalize input stream
    if in_fasta == '-':
        _input = sys.stdin
    else:
        _input = open(in_fasta)

    fasta_out = FastaIO.FastaWriter(sys.stdout, wrap=None)
    fasta_out.write_file(parse(_input, 'fasta'))
