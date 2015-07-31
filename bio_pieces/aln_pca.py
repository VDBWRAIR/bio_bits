from __future__ import print_function

import itertools
from collections import OrderedDict

from Bio import SeqIO
import pandas as pd
import numpy as np

def pairwise_identity(seq1, seq2, invalid_chars=None):
    '''
    Compare two sequences by counting only bases that
    are identical between them

    :param str seq1: string of chars
    :param str seq2: strin of chars
    :return: sum of same base positions
    :raises: ValueError if any sequence contains any invalid_chars or
             if sequence lengths are not identical
    '''
    if len(seq1) != len(seq2):
        raise ValueError('Sequence lengths did not match')

    def invalid(c):
        if invalid_chars is None:
            return
        if c in invalid_chars:
            print("c: {0} invalid_chars: {1}".format(c,invalid_chars))
            raise ValueError('{0} is an invalid character'.format(c))

    ident = 0
    for x,y in itertools.izip(seq1,seq2):
        invalid(x)
        invalid(y)
        if x.lower() == y.lower():
            ident += 1
    return ident

def index_fasta(aln_fh):
    '''
    Return a pandas.Series for the fasta sequence alignment in order
    of the sequences in the file

    :param file aln_fh: fasta file like iterator
    :return: pandas.Series indexed by Bio.SeqRecord.id
    '''
    # Build list of tuples (id,seq) and convert to series later
    seq_index = OrderedDict()
    for record in SeqIO.parse(aln_fh, 'fasta'):
        seq_index[record.id] =  str(record.seq)
    return pd.Series(seq_index)

def identity_matrix(aln):
    '''
    Build an identity matrix from all the pairwise identities
    of all sequences in the supplied fasta index

    n^2 loop over sequences to generate all identities.
    Skips identies of sequences against themselves
    Only does top right calculations of matrix and copies values
    into bottom left as they are identical

    :param mapping aln: Aligned indexed fasta
    :return: pandas.DataFrame representing identity matrix
    '''
    id_matrix = np.empty([len(aln),len(aln)])
    for i in range(len(aln)):
        for j in range(len(aln)):
            #print(i,j)
            # Don't need to compute identity against itself
            if i == j:
                #print("Using length")
                id_matrix[i][j] = len(aln[j])
            # Only compute top right of matrix
            elif i > j:
                #print("Copying {0}{1} from {1}{0}".format(j,i))
                id_matrix[i][j] = id_matrix[j][i]
            else:
                _id = pairwise_identity(aln[i], aln[j])
                #print("Pident of {0} and {1} is {2}".format(aln[i],aln[j],_id))
                id_matrix[i][j] = _id
    return pd.DataFrame(id_matrix, index=aln.keys(), columns=aln.keys())
