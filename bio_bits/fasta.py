from __future__ import print_function
import sys
import argparse
import string
import functools

from Bio.SeqIO import parse, FastaIO
from Bio.Data import IUPACData

from . import util
from . import ctleptop 

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
    parser.add_argument(
        '--disambiguate', action='store_true', default=False,
        help='Generate all permutations sequence based on ambiguous nucleotide'
            ' sequences'
    )
    return parser.parse_args()

def writefasta(records, outfile=sys.stdout, wrap=None):
    '''
    :param SeqRecord records: Sequence of SeqRecord
    :param handle|string outfile: Output file or handle to write to
    :param bool wrap: Wrap sequences at 80 chars or not
    '''
    fasta_out = FastaIO.FastaWriter(outfile, wrap=wrap)
    fasta_out.write_file(records)

def splitfasta(records, wrap=False):
    for seq in records:
        outfile = str(seq.id)
        for p in string.punctuation:
            outfile = outfile.replace(p, '_')
        outfile = outfile + '.fasta'
        with open(outfile, 'w') as fh:
            writefasta([seq], fh, wrap)

def disambiguate(records):
    all_records = []
    for record in records:
        make_rec = functools.partial(
            util.make_seqrec,
            id=record.id, name=record.name, description=record.description
        )
        # Get the input sequence
        sequence = str(record.seq)
        # Generate all permutations of that sequence
        #perms = list(ctleptop.nearbyPermutations(sequence))
        sys.stderr.write(record.id + '\n')
        perms = permutate_ambiguous_sequence(sequence)
        for perm in perms:
            yield make_rec(perm)

def permutate_ambiguous_sequence(seq_str):
    '''
    '''
    # Abiguous mapping table from BioPython
    amb_values = IUPACData.ambiguous_dna_values
    a_bases = []
    total_perms = 1
    for nt in seq_str:
        amb_bases = amb_values.get(nt, nt)
        if len(amb_bases) > 1:
            a_bases.append(nt)
            total_perms *= len(amb_bases)
    sys.stderr.write("Sequence has {0} ambiguous bases over {1} total bases for a total of {2} permutations\n".format(len(a_bases), len(seq_str), total_perms))
    # Start ambiguous sequences with our input sequence
    amb_seqs = [seq_str]
    # i holds current position
    for i in range(len(seq_str)):
        nt = seq_str[i]
        amb_bases = amb_values.get(nt, nt)
        # Skip all non-ambiguous bases
        if len(amb_values) == 1:
            continue
        #print("i: {0}".format(i))
        cur_seqs = []
        # Go through each sequence again and and generate ambiguous bases
        for seq in amb_seqs:
            #print("nt: {0}".format(nt))
            # build up permutations for the current ambiguous base
            for base in amb_bases:
                #print("base: {0}".format(base))
                #print("seq[:i] + base + seq[i+1:]".format(seq[:i], base, seq[i+1:]))
                cur_seqs.append(seq[:i] + base + seq[i+1:])
                #print("cur_seqs: {0}".format(cur_seqs))
        amb_seqs = cur_seqs
        #print("amb_seqs: {0}".format(amb_seqs))
    return amb_seqs

def main():
    args = parse_args()
    in_fasta = args.fasta

    # Normalize input stream
    if in_fasta == '-':
        _input = sys.stdin
    else:
        _input = open(in_fasta)
    input_records = parse(_input, 'fasta')

    if args.split:
        splitfasta(input_records, args.wrap)
    elif args.disambiguate:
        recs = disambiguate(input_records)
        writefasta(recs, sys.stdout, args.wrap)
    else:
        writefasta(input_records, sys.stdout, wrap=args.wrap)
