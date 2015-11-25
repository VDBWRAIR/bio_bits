from functools import partial
from itertools import dropwhile, takewhile, starmap, combinations, imap, ifilter
import itertools
from toolz.itertoolz import  accumulate, first, last, drop, mapcat, take
from toolz.functoolz import complement, compose
from operator import itemgetter
import fn
import sh
import re

#from operator import xor
#def ilen(seq): sum(1 for _ in seq)
#indel_pos = compose(ilen, partial(takewhile, '-'.__ne__))
#def first_indel_pos(ref, query, cigar_string):
#    ralign, qalign = undo_cigar(ref, query, cigar_string)
#    ins, _del = indel_pos(ralign), indel_pos(qalign)
#    assert xor(ins, _del), "Indel contains insertion and deletion"
#    return ins if ins else _del
#reduce(next_cigar_position, cig_elements, (0,0))
# read = lambda x: vim.current.buffer.append(str(x))

plus_n = lambda x, n : x + n
zero = lambda x, y : x
# changing plus_n to alter strings won't quite work for the string transformations because
# adding gaps is dependent on there being a difference between ref and pos.
dx_dy={
    'M' : (plus_n, plus_n), 'X' : (plus_n, plus_n), '=' : (plus_n, plus_n),
    'D' : (plus_n, zero), 'plus_n' : (plus_n, zero),
    'I' : (zero, plus_n), 'S' : (zero, plus_n),
    'H' : (zero, zero)
}

def transform(seqs, positions):
    (ref, query), ((rpos, qpos), (a, b)) = seqs, positions
    def fill(s, pos): return s[:pos] + '-'*abs(a-b) + s[pos:]
    if a < b: ref = fill(ref, rpos)
    elif a > b: query = fill(query, qpos)
    return  ref, query

def next_cigar_position(positions, cigarval, withtype=False):
    if withtype: (count, sym), (refpos, qpos, oldcount, oldsym) = cigarval, positions
    else: (count, sym), (refpos, qpos) = cigarval, positions
    f = dx_dy[sym]
    if withtype: return f[0](refpos, count), f[1](qpos, count), count, sym
    return f[0](refpos, count), f[1](qpos, count)
cig_diffs = partial(map, partial(next_cigar_position, (0, 0)))

def undo_cigar(x, y, z):
    cig = get_cig(z) #return reduce(transform, zip(gen_cigar(z), cig_diffs(z)), (x, y))
    return reduce(transform, zip(_gen_cigar(cig), cig_diffs(cig)), (x, y))

#in order to *totally* recreate the reads, you need to account for the other weird tag


split_cig = re.compile(r'(?:([0-9]+)([A-Z]))').findall
get_cig = lambda C: map(lambda x: (int(x[0]), x[1]), split_cig(C))
_gen_cigar = lambda c: accumulate(next_cigar_position, c, (0, 0))
gen_cigar = compose(_gen_cigar, get_cig)


def get_insert_quals(ref, pos, bases):

    '''separate i/o and destructuring from logic.
    make generic for indel/snp.'''
    view = sh.samtools('view', "{ref}:{pos}-{pos}")
    #split tabs, etc.
    has_insert = lambda l: 'I' in l[5]
    # etc.
    #samtools calcmd?
# ref first
assert undo_cigar('AAGC', 'AGTT', '1M1D1M1X1I') == ('AAGC-', 'A-GTT')
assert list(starmap(undo_cigar, [
["AGG", "AG"   ,"1M1D1M"],
["AGG", "AG"   ,"1M1D1M"],
["TT" , "TAT"  , "1M1I1M"],
["TA" , "TATAA", "1M3I1M"],
['AAGC', 'AGTT', '1M1D1M1X1I']])) ==  [('AGG', 'A-G'), ('AGG', 'A-G'), ('T-T', 'TAT'), ('T---A', 'TATAA'), ('AAGC-', 'A-GTT')]
assert list(gen_cigar("1M2D1M")) == [(0, 0), (1, 1), (3, 1), (4, 2)]



def intersection(pred, seq): return last(takewhile(complement(pred), seq))
def firstwhere(pred, seq): return first(dropwhile(complement(pred), seq))

assert intersection(lambda x: x[0] > 2, gen_cigar("1M2D1M")) == (1, 1)
assert firstwhere(lambda x: x[0] >= 2, gen_cigar("1M1D1M")) ==  (2, 1)

#firstwhere(lambda x: x[1] == 'D', gen_cigar("1M1D1M"))
#firstwhere(lambda x: x[1] == 'D', get_cig("1M1D1M"))
c = "1M1D1M"
''' corresponds to ('AAGC-', 'A-GTT')'''
assert reduce(next_cigar_position,  \
       takewhile(lambda x: x[1] != 'I', get_cig('1M1D1M1X1I')), (0, 0)) == (4, 3)

def mutpos(cig_str, types):
    ''':return (refpos, mutpos'''
    assert any(map(types.__contains__, cig_str)), "Cigar string %s does not have mutation of type %s" % (cig_str, types)
    return reduce(next_cigar_position, takewhile(lambda x: x[1] not in types, get_cig(cig_str)), (0, 0))

foo = lambda c, t: reduce(next_cigar_position, filter(lambda x: x[1] in t, get_cig(c)), (0, 0))
(lambda c, t: accumulate(next_cigar_position, filter(lambda x: x[1] in t, get_cig(c)), (0, 0)))
indelpos = partial(mutpos, types='DI')
query_indel_pos = compose(last, indelpos)
ref_indel_pos = compose(first,  indelpos)

ref_del_pos = compose(first, partial(mutpos,  types='D'))
ref_ins_pos = compose(first, partial(mutpos,  types='I'))
ref_snp_pos = compose(first, partial(mutpos,  types='X'))


def mutations(cig_str):
    return drop(1, accumulate(partial(next_cigar_position, withtype=True), get_cig(cig_str), (0, 0, 0, 0)))

snps = compose(partial(filter, lambda x: x[-1] == 'X'), mutations)

def snp_rpos_base(qstring, cig_str):
    return compose(partial(mapcat, partial(rpos_base, qstring)), snps)(cig_str)

def rpos_base(qstring, (rpos, qpos, count, type)):
    ''':return refpos, snp'''
    return ((qstring[(qpos-count)+i], rpos +i) for i in xrange(0, count))


#wishful thinking about complex
#{'ins' : 'I', 'del' : 'D', 'snp' 'X', 'complex' : 'X'}

assert ref_ins_pos('1M1D1M1X1I') == 4
#use assertRaisesRegex
try:
    mutpos('1M1D1M', 'X')
except AssertionError:
    pass

#    cigs = get_cig(CIGAR)
#    indel_effect = sum(map(first, filter(compose('D'.__eq__, last), cigs))) - sum(map(first, filter(compose('I'.__eq__, last), cigs)))

def get_info(samrow):
    #extract query_string and cigar_string from file
    # have to add map positions to computed rpos's.
    qname, mpos, cigar, qstring = itemgetter(0, 3, 5, 9)(samrow.split('\t'))
    mpos = int(mpos)
    _snps = snp_rpos_base(qstring, cigar)
    return map(lambda (b,i):(b, i+mpos), _snps)

# split by RG, then get_info, want to keep track of qname to prob.
#raw_longs, raw_shorts = fn.iters.partition(
lhn = ghn = 3
#shareCount = compose(len, set.intersection)
shareCount = compose(len, lambda (x,y): x & y)
is_hap = lambda n: lambda x: shareCount(x) > n
# have to get cigar string with 'X' in it.
if __name__ == '__main__':
    import multiprocessing  as mp
    import sys
    _fn = sys.argv[1]
    p = mp.Pool()
    shortReads = imap(get_info, take(3, sh.samtools.view(_fn, _iter=True)))
    shortReads = imap(set, shortReads)
    #import ipdb; ipdb.set_trace();
    cartProduct = itertools.product(*list(itertools.tee(shortReads)))
    #print map(lambda (x,y): (type(x), type(y)), cartProduct)
    shortReads = ifilter(is_hap(lhn), cartProduct)
    for r in shortReads:
        if r: print r


#ghs = filter(is_hap(ghn), combinations(longReads, longReads))
#lhs = filter(is_hap(lhn), combinations(shortReads, shortReads))

