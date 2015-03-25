# -*- coding: utf-8 -*-

import re
from Bio import SeqIO
from collections import defaultdict

def make_genemap( sequencefile, splitter, colnum, format='fasta' ):
    genemap = defaultdict(list)
    for record in SeqIO.parse(sequencefile, format):
        parts = re.split(splitter,record.id)
        gene = parts[colnum]
        genemap[gene].append( record )

    return genemap

def write_split_files( genemap, outformat='fasta' ):
    for gene, records in genemap.iteritems():
        SeqIO.write(records, open(gene+'.fasta','w'), outformat)

def parse_args():
    from argparse import ArgumentParser
    
    parser = ArgumentParser(
        description='Splits a fasta file into multiple files by ' \
            'grouping sequences by a specified column in the sequence ' \
            'identifier'
    )

    parser.add_argument(
        'seqfile',
        help='Sequence file path to split'
    )

    parser.add_argument(
        '--delimiter',
        '-d',
        default='__',
        help='The delimiter in the sequence id line to use'\
            '[Default: %(default)s]'
    )

    parser.add_argument(
        '--colnum',
        '-c',
        default=2,
        type=int,
        help='The column in the sequence identifier line to use to group'\
            '[Default: %(default)s]'
    )

    parser.add_argument(
        '--file-type',
        '-ift',
        dest='seqfiletype',
        default='fasta',
        help='The file type the input file is[Default: %(default)s]'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    genemap = make_genemap(
        args.seqfile,
        args.delimiter,
        args.colnum,
        args.seqfiletype
    )
    write_split_files( genemap, args.seqfiletype )
