try:
    import unittest2 as unittest
except ImportError:
    import unittest

import mock
import sys
from bio_pieces import amos, myargs
from test_amos import make_red, make_tle, make_ctg, make_amos_string
import os
#TODO: Restore a second input file.
#TODO: Figure out a better way to produce test inputs. SeqIO? 
'''
Assumptions: 
    1. Each RED corresponds to exactly one FASTQ entry.
Want to test that we handle cases like:
    1. The same read maps to two contigs
    2. Two REDS/fastq entries have the same sequence string

'''
THISD= '.'
inputdir = os.path.join(THISD, 'testinput')

class TestAmos2Fastq(unittest.TestCase): 
#           def setUp(self):
#               self.redlist, self.tlelist, self.ctglist, self.amosstr = make_amos_string()
#               #with open('tests/foo.afg', 'w') as outfile: outfile.write(self.amosstr)
#               self.hashed={'ATCG' : [RED(i, i, 'ATCG', 'D'*4) for i in (1, 2, 9)],
#                       'TCAG' : [RED(3, 3, 'TCAG', 'D'*4)], 
#                       'ACTG' : [RED(4, 4, 'ACTG', 'D'*4)],
#                       'GTCA' : [RED(5, 5, 'GTCA', 'D'*4)],
#                       'GACT' : [RED(6, 6, 'GACT', 'D'*4)],
#                       'GTAC' : [RED(7, 7, 'GTAC', 'D'*4)],
#                       'CTGA' : [RED(8, 8, 'CTGA', 'D'*4)]
#                           }
#               fastqs = ['foo1.fastq', 'foo2.fastq']
#               self.fastqs = [os.join(inputdir, fq) for fq in fastqs]
#               pass

#    def test_hash_reds_by_sequence(self):
#        result = a2f.hash_reds_by_sequence(self.redlist)
#        self.assertEquals(self.hashed, result) 
#        pass
#
#    def test_iter_multi_fastq(self):
#        result = a2f.multi_fastq_iterator(self.fastqs)
#        self.assertEquals(10, len(result))
#
#    def test_get_fastq_as_df(self):
#        result_df =a2f.fastq_as_df( open('expected/expectedfoo1.fastq'))
#
#
#        pass
    def test_main(self):
        sys.argv='../bio_pieces/myargs.py testinput/foo1.fastq testinput/foo2.fastq --amos testinput/foo.afg'.split(' ')
        myargs.main()


