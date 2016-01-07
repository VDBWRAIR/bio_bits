from __future__ import print_function
import sys
import argparse

from Bio.SeqIO import parse

def parse_args():
    parser = argparse.ArgumentParser(
        '''Simplistic fasta manipulator'''
    )
    parser.add_argument(
        'fasta', default='-',
        help='Fasta file path or - to read from standard input'
    )
    return parser.parse_args()

def main():
    in_fasta = sys.argv[1]

    if in_fasta == '-':
        _input = sys.stdin
    else:
        _input = open(in_fasta)

    for rec in parse(_input, 'fasta'):
        sys.stdout.write('>{0}\n'.format(rec.description))
        sys.stdout.write(str(rec.seq) + '\n')
        sys.stdout.flush()
