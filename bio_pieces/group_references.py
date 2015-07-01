'''
Usage: group_references <samfile> --outdir <DIR>

Options:
    --outdir=<DIR>,-o=<DIR>   output directory [Default: parse_contigs_out]

Create separate fastq files for each reference in a samfile.
'''

from docopt import docopt
import pandas as pd
from schema import Schema, Use
import sys
import os
import sh
if sys.version[0] == '3':
    from io import StringIO as BytesIO
else:
    from io import BytesIO
#samview = check_output('samtools view {0}'.format(samfile))

sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
def samview_to_df(rawtext):
    samtext = '\n'.join( fixline(row) for row in rawtext.split('\n') )
    as_bytes = BytesIO(samtext)
    #Write test to ensure handles varaibles #columns
    return pd.read_csv(as_bytes, names=sam_columns, usecols=sam_columns, delimiter='\t', header=None, squeeze=True)

def fixline(row):
    newrow = []
    cols = row.split('\t')[:len(sam_columns)]
    if len(cols) > 1 and cols[2] == '*':
        cols[2] = 'unmapped'
    return '\t'.join(cols)

def get_seqs_by_ctg(outdir, rawtext):
    sam_df = samview_to_df(rawtext)
    contig_groups = sam_df.groupby('RNAME')
    fastq = "@{0}\n{1}\n+\n{2}".format
    for group in contig_groups:
        ref, reads = group[0], group[1]
        with open("{0}/{1}.group.fq".format(outdir, ref), 'w') as out:
            map(out.writelines, '\n'.join(map(fastq, reads.QNAME, reads.SEQ, reads.QUAL)))
            out.write('\n')


def main():
    raw_args = docopt(__doc__)
    scheme = Schema({
        '<samfile>' : str,
        '--outdir' : str})
    parsed_args = scheme.validate(raw_args)
    outdir = parsed_args['--outdir']
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sam = str(sh.samtools('view', parsed_args['<samfile>']))
    get_seqs_by_ctg(outdir, sam)
    return 0

if __name__ == '__main__':
    main()
