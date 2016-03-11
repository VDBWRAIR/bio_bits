from . import unittest

from bio_bits import util

class TestMakeSeqrec(unittest.TestCase):
    def test_makes_record(self):
        r = util.make_seqrec('ATGC', 'id', 'name', 'description')
        self.assertEqual(str(r.seq), 'ATGC')
        self.assertEqual(r.id, 'id')
