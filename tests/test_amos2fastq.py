try:
    import unittest2 as unittest
except ImportError:
    import unittest

import mock
import sys
from bio_pieces import amos, myargs, amos2fastq_main
from test_amos import make_red, make_tle, make_ctg, make_amos_string
import os
from functools import partial
#TODO: Restore a second input file.
#TODO: Figure out a better way to produce test inputs. SeqIO?
'''
Assumptions:
    1. Each RED corresponds to exactly one FASTQ entry.
Want to test that we handle cases like:
    1. The same read maps to two contigs
    2. Two REDS/fastq entries have the same sequence string

'''
THISD = os.path.dirname(os.path.abspath(__file__))
inputdir = os.path.join(THISD, 'testinput')
outputdir = os.path.join(THISD, '.')
expecteddir  = os.path.join(THISD, 'expected')

class TestAmos2Fastq(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_seqs_by_ctg(self):
        pass

    def test_extract_dfs_by_iids(self):
        pass


class TestAmos2FastqMain(unittest.TestCase):

    def setUp(self):
        input_join = partial(os.path.join, inputdir)
        output_join = partial(os.path.join, outputdir)
        expected_join = partial(os.path.join, expecteddir)
        outnames = ['contig-0.fastq', 'contig-1000000.fastq', 'contig-2000000.fastq']
        self.fastqs = map(input_join, ['1.R1.unmap.fastq', '1.R2.unmap.fastq'])
        self.amos = input_join('regression.afg')
        self.expected_files = map(expected_join, outnames)
        self.result_files = outnames #map(output_join, outnames)

    @mock.patch('bio_pieces.amos2fastq_main.docopt')
    def test_main(self, mockopt):
        mockopt.return_value = {'<fastqs>' : self.fastqs, '--amos' : self.amos}
        read_file = lambda fn: open(fn).read()
        amos2fastq_main.main()
        map(self.assertTrue, map(os.path.exists, self.result_files))
        map(self.assertEquals, *map(read_file, self.expected_files, self.result_files))
    pass

#    def test_main(self):
#        sys.argv='../bio_pieces/myargs.py testinput/foo1.fastq testinput/foo2.fastq --amos testinput/foo.afg'.split(' ')
#        myargs.main()

class TestPandasUtils(unittest.TestCase):

    def setUp(self):
        pass



