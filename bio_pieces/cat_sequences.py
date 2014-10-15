import re
from os.path import basename
import sys

import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main():
    args = parse_args()
    srpf = sort_seq_files(args.seqfiles, args.delimiter, args.sortkeys)
    records = combine_seqs_inorder(
        srpf,
		args.delimiter,
		args.sortkeys,
    )
    SeqIO.write(records, sys.stdout, 'fasta')

def parse_args():
    parser = argparse.ArgumentParser(
        description = 'Merges multiple fasta files into one list of '\
            'concatenated single records.' \
            'Essentially creates a 2-D matrix from all the files and '\
            'their records. Then sorts all columns(files). Then sorts '\
            'all rows(sequence records). Once sorting is done, then '\
            'concatenates across the rows in ascending order. This '\
            'is useful for building genomes from segment fasta files.'
    )

    parser.add_argument(
        'seqfiles',
        nargs='+',
        help='List of sequence files to concatenate'
    )

    parser.add_argument(
        '--delimiter',
        '-d',
        default='__',
        help='The delimiter to use to split the sequence id[Default: %(default)s]'
    )

    parser.add_argument(
        '--sortkeys',
        '-k',
        default=[1,2],
        nargs='+',
        help='Which columns in the sequence id to use after split using ' \
            'delmiter[Default: %(default)s]'
    )

    return parser.parse_args()

def sort_seq_files(seqfiles, delimiter, sortkeys):
    '''
    Create a dictionary mapping each seqfilename to its sequences
    Sort the sequences in each item of the dictionary

    Return dictionary with sorted seqrecords for each file
    '''
    seqsperfile = {}
    for seqfile in seqfiles:
        seqs = list(SeqIO.parse(seqfile, 'fasta'))
        sort_sequences(seqs, delimiter, sortkeys)
        seqsperfile[basename(seqfile)] = seqs
    return seqsperfile

def combine_seqs_inorder(seqrecordsperfile, delimiter, sortkeys):
    '''
    Assumes seqrecordsperfile has each seqrecord list already sorted
    Can use sort_seq_files to get this dictionary

    Go through each item in seqrecordsperfile and sort all first
    elements, then sort all second elements
    When finished sorting an a list of first elements then join them
    in the sorted order

    seqrecordsperfile - dictionary mapping {name: [seqrecords],}
        probably from sort_seq_files
    delmiter - What to split each seqrecord.id on
    sortkeys - What columns after splitting to sort on in order

    return single seqrecord list with combined seqrecords
    '''
    # This will hold the resulting concatted sequence records
    seqrecords = []
    # Gets only the seqrecords lists which should be sorted now
    perfilerecords = seqrecordsperfile.values()
    # Ensure all list lengths are same
    if len(set(map(len,perfilerecords))) != 1:
        # We will handle this better some time later
        raise Exception("Not all sequence record lists have same length")
    # Transpose lists such that we get lists per index across them all
    transposed = map(list, zip(*perfilerecords))
    # Sort each of the transposed lists in place
    # then cat them to get our resulting seqrecords
    for i in range(len(transposed)):
        sort_sequences(transposed[i], delimiter, sortkeys)
        rec = cat_seqrecords(transposed[i], delimiter, sortkeys)
        seqrecords.append(rec)
    return seqrecords

def cat_seqrecords(seqrecords, delimiter, sortkeys, keepdescriptions=True):
    '''
    Simply append from left to right the seqrecord.seq.seq in each element
    and return the resulting Bio.seqrecord
    The id of the seqrecord will be the first item of sortkeys
        as that should be the same hopefully in each seqrecord
    The description will be all id's joined together by comma unless keepdescriptions
        is False then it will be ''
    '''
    newseq = []
    newid = ''
    newdescription = []
    for rec in seqrecords:
        splitid = split_seq_id(rec, delimiter)
        # Again sortkeys is 1-indexed
        newid = splitid[sortkeys[0]-1]
        newseq.append(str(rec.seq))
        if keepdescriptions:
            newdescription.append(rec.id)
    # Set up newdescription
    if not newdescription:
        newdescription = ''
    else:
        newdescription = ','.join(newdescription)
    # New sequence record
    # Has same alphabet as first item in original list
    # Has concatted sequences
    newrecord = SeqRecord(
        Seq(
            ''.join(newseq),
            seqrecords[0].seq.alphabet
        ),
        id = newid,
        description = newdescription,
    )

    return newrecord

def sort_sequences(seqs, delimiter, sortkeys):
    '''
    Sort a list of Biopython seqrecords using sortkeys after splitting each
    seq.id with the delimiter

    Items are sorted in place
    '''
    def keyfunc(seqrec):
        '''Join sortkeys into a string'''
        ids = split_seq_id(seqrec, delimiter)
        key = ''
        for k in sortkeys:
            # sortkeys is 1-indexed not 0-indexed
            key += ids[k-1]
        return key
    seqs.sort(key=keyfunc)

def split_seq_id(seqrecord, delimiter):
    '''
    Simply split seqrecord.id using delimiter and return the list

    returns list of split items
    '''
    return re.split(delimiter, seqrecord.id)
