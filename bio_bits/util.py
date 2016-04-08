from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet

def make_seqrec(*args, **kwargs):
    '''
    Just call SeqRecord constructor, but wraps Seq object around string
    for you

    make_seqrec('ATGC', 'id', 'name', 'description')
    '''
    return SeqRecord(Seq(args[0]), *args[1:], **kwargs)
