import sys
from . import unittest

from bio_bits import fasta, util
import mock

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

class TestPermutateAmbiguousSequence(unittest.TestCase):
    def test_works_no_ambig(self):
        r = fasta.permutate_ambiguous_sequence('foo', 'ATG-C')
        self.assertEqual(['ATG-C'], r)

    def test_works_all_ambig(self):
        r = fasta.permutate_ambiguous_sequence('foo', 'RRR')
        self.assertEqual(
            sorted(['AAA', 'AAG', 'AGA', 'GAA', 'AGG', 'GGA', 'GGG', 'GAG']), sorted(r)
        )

    def test_known_test(self):
        r = fasta.permutate_ambiguous_sequence('foo', 'ATGRB')
        self.assertEqual(
            sorted(['ATGAT', 'ATGAG', 'ATGAC', 'ATGGT', 'ATGGG', 'ATGGC']), sorted(r)
        )

    def test_over_max_perms(self):
        with mock.patch.object(sys, 'stderr') as m_stderr:
            # Should generate 128 permutations
            r = fasta.permutate_ambiguous_sequence('foo', 'R'*7)
            m_stderr.write.assert_called_with('Sequence foo has 7 ambiguous bases that would produce 128 permutations and was skipped\n')
            self.assertEqual([], r)
