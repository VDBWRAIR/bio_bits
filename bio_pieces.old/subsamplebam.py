import sys
import argparse
import random

import pysam

class MinimumDepthException(Exception):
    ''' Exception for when samplesize > depth '''

def get_readset_for_reference(samfile, referenceid, samplesize):
    '''
    For a given reference id, fetch the unique set of read names that guarantees
    samplesize depth across the reference

    :raises MinimumDepthException: if any position on the reference does not contain
    enough reads to fill samplesize depth

    :param pysam.AlignmentFile samfile: samfile to get pileup to iterate
    :param str referenceid: reference name to use for pileup call
    :param int samplesize: subsample depth

    :return: Unique set of reads that will create samplesize depth across reference
    :rtype set:
    '''
    # Read ids to keep
    keepset = set([])
    # Read ids that are have been seen at a position so
    # no longer can be used
    excludeset = set([])

    # Iterate left->right across the reference
    for pilecol in samfile.pileup(reference=referenceid):
        # Will contain set of all reads at this position
        readset = set([])
        # get all read names at this position
        for read in pilecol.pileups:
            readset.add(read.alignment.query_name)
        # items that are already in keepset
        commonnames = keepset.intersection(readset)
        print 'Common: {0}'.format(commonnames)
        # only items that can be selected as new
        leftover = readset - commonnames - excludeset
        print 'leftover: {0}'.format(leftover)
        # How many reads we still need to add for this position
        randreads = samplesize - len(commonnames)
        print 'randreads: {0}'.format(randreads)
        # randomly selected read names
        #selected = random.sample(leftover, randreads)
        selected = list(leftover)[:samplesize]
        print 'selected: {0}'.format(selected)
        # Add to the master list
        keepset.update(selected)
        # We can't add the others later on
        excludeset.update(leftover - selected)
    return keepset

def write_subsampled_sam(input, output, samplesize):
    '''
    Write subsampled bam reads as SAM format to outputfh
    Each position on each reference will contain exactly samplesize
    depth
    '''
    inputsam = pysam.AlignmentFile(input, 'rb')
    for refname in inputsam.references:
        usereads = get_readset_for_reference(inputsam, refname, samplesize)

def main():
    args = parse_args()
    write_subsampled_sam(args.input, args.soutput, args.samplesize)

def parse_args():
    parser = argparse.ArgumentParser(''' 
        Sub-samples a bam file such that there is even coverage across the genome by selecting
        random reads that create even coverage.
    ''')

    parser.add_argument(
        'input',
        help='Input bam file path or - for stdin'
    )

    parser.add_argument(
        'output',
        default=sys.stdout,
        help='Sub-sampled bam output as SAM format. Either file path or - for stdout'
    )

    parser.add_argument(
        'samplesize',
        help='The wanted coverage depth'
    )

    return parser.parse_args()
