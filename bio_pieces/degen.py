from __future__ import print_function
from Bio import Entrez, SeqIO, GenBank
from functools import partial
from collections import namedtuple
from itertools import starmap, product, imap as map
import operator
import re
Entrez.email = "micheal.panciera.work@gmail.com"

def compose2(f, g):
    def inner(*args, **kwargs):
        return f(g(*args, **kwargs))
    return inner
def compose(*funcs): return reduce(compose2, funcs)

DEGENS = ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M']
Gene = namedtuple('Gene', ['name', 'start', 'stop'])


#TODO: but `genes` includes CDS which is whole genome. I guess that's not the worst thing though
def unfinished():
    rec = fetch_record_by_id(id)
    res = GenBank.RecordParser().parse(StringIO.StringIO(rec))
    codes = filter(lambda x: x.key == 'mat_peptide' or x.key == 'CDS', res.features)
    genes = filter( lambda x: filter(lambda y: y.key == '/product=', x.qualifiers), codes)

def gene_from_line(s): return Gene(get_name(s), *get_interval(s) )
def fetch_record_by_id(_id): return Entrez.efetch(db="nucleotide", id=_id, rettype='gb', retmode='text').read()

gene_lines = re.compile(r'(mat_peptide[^\n]+\n[^\n]+)').findall
raw_nums = lambda x: re.compile(r'([0-9]+)\.\.([0-9]+)\n').search(x).groups()
get_interval = compose(partial(map, int), raw_nums)
get_name = lambda x: re.compile(r'/product="([^"]+)').search(x).groups()[0]
degen_positions = lambda seq:  (m.start() for m in re.finditer('|'.join(DEGENS), seq))

genes_form_text = compose(partial(map, gene_from_line), gene_lines)
fetch_genes_by_id = compose(genes_form_text, fetch_record_by_id)

rowformat='{0}\t{1}\t{2}'

def get_gene_pos_seq(genes, seq):
    ''' :return list of tuples of form: (gene name, position, nt)... '''
    _degen_positions = degen_positions(str(seq))
    perms = product(genes, _degen_positions)
    result = ((gene.name, pos, seq[pos]) for gene, pos in perms if  gene.start <= pos <= gene.stop)
    return result

def get_result_info(infasta, ref_id=None, gene_file=None):
    assert operator.xor(bool(ref_id), bool(gene_file)), "Must specify reference id or gene annotation file, but not both. \
        You specified reference id %s and gene file %s" % (ref_id, gene_file)
    if ref_id:
        genes = fetch_genes_by_id(ref_id)
    else:
        raise NotImplementedError('Currently only works with GenBank reference ids.')
    fasta = SeqIO.parse(infasta, format='fasta')
    return map(partial(get_gene_pos_seq, genes), fasta)

pretty_table = compose('\n'.join, partial(starmap, rowformat.format))
get_result_table = compose(partial(map, pretty_table), get_result_info)
