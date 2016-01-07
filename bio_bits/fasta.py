import sys

from Bio.SeqIO import parse

in_fasta = sys.argv[1]

if in_fasta == '-':
    _input = sys.stdin
else:
    _input = open(in_fasta)

for rec in parse(_input, 'fasta'):
    sys.stdout.write('>{0} {1}\n'.format(rec.id, rec.description))
    sys.stdout.write(str(rec.seq) + '\n')
    sys.stdout.flush()
