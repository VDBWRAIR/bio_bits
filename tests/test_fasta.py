from . import unittest

from bio_bits import fasta, util

class TestDisambiguate(unittest.TestCase):
    def sort_seqrecs(self, recs):
        return sorted(recs, key=lambda rec: rec.seq)

    def assertSeqRecEqual(self, rec1, rec2):
        attrs = ['seq', 'id', 'name', 'description']
        for attr in attrs:
            self.assertEqual(
                getattr(rec1, attr), getattr(rec2, attr)
            )

    def assertSeqRecsEqual(self, seq1, seq2):
        for s1, s2 in zip(seq1, seq2):
            self.assertSeqRecEqual(s1, s2)

    def test_builds_records(self):
        start_seq = util.make_seqrec('ATGRB')
        all_seqs = ['ATGAC', 'ATGAG', 'ATGAT', 'ATGGC', 'ATGGG', 'ATGGT']
        all_seqrecs = map(util.make_seqrec, all_seqs)
        r = fasta.disambiguate([start_seq])
        self.assertSeqRecsEqual(self.sort_seqrecs(all_seqrecs), self.sort_seqrecs(r))
