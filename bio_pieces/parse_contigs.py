'''
Usage: group_references <samfile>

Create separate fastq files for each reference in a samfile.
'''

from docopt import docopt
import pandas as pd
import io
from schema import Schema, Use

#samview = check_output('samtools view {0}'.format(samfile))

sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
def samview_to_df(samtext):
    as_bytes = io.BytesIO(samtext)
    return pd.read_csv(as_bytes, names=sam_columns, usecols=sam_columns, delimiter='\t')

def get_seqs_by_ctg(samtext):
    sam_df = samview_to_df(samtext)
    contig_groups = sam_df.groupby('RNAME')
    fastq = "{0}\n{1}\n+\n{2}".format
    for group in contig_groups:
        ref, reads = group[0], group[1]
        with open("{0}.group.fq".format(ref), 'w') as out:
            map(out.writelines, '\n'.join(map(fastq, reads.QNAME, reads.SEQ, reads.QUAL)))

def main():
    raw_args = docopt(__doc__)
    scheme = Schema({'<samfile>' : Use(open, error='Samfile must be readable')})
    parsed_args = scheme.validate(raw_args)
    get_seqs_by_ctg(parsed_args['<samfile>'].read())
    return 0

if __name__ == '__main__':
    main()

