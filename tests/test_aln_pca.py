from StringIO import StringIO

from . import unittest, THIS
import pandas as pd

from bio_pieces import aln_pca

aln = [
    ('>id1','ATGC'),
    ('>id2','AAGC'),
    ('>id3','CGTR'),
    ('>id4','ACGT')
]

id_matrix = pd.DataFrame(
    {
        'id1': [4.,3.,0.,2.],
        'id2': [3.,4.,0.,2.],
        'id3': [0.,0.,4.,0.],
        'id4': [2.,2.,0.,4.]
    },
    index=['id1','id2','id3','id4']
)

class TestPairwiseIdentity(unittest.TestCase):
    def test_id_same(self):
        seq = 'A' * 10
        r = aln_pca.pairwise_identity(seq, seq)
        self.assertEqual(len(seq), r)

    def test_id_complete_diff(self):
        seq = 'A' * 10
        r = aln_pca.pairwise_identity(seq, seq.replace('A','C'))
        self.assertEqual(0, r)

    def test_id_half_diff(self):
        r = aln_pca.pairwise_identity('A'*10, 'A'*5+'C'*5)
        self.assertEqual(5, r)

    def test_diff_sizes(self):
        self.assertRaises(
            ValueError,
            aln_pca.pairwise_identity, 'A'*5, 'A'*4
        )

    def test_ignores_case(self):
        r = aln_pca.pairwise_identity('a'*10, 'A'*10)
        self.assertEqual(10, r)

    def test_non_seq_char_dash(self):
        self.assertRaises(
            ValueError,
            aln_pca.pairwise_identity, 'A-CG', 'AACG', '- '
        )

    def test_non_seq_char_space(self):
        self.assertRaises(
            ValueError,
            aln_pca.pairwise_identity, 'AACG', 'A CG', '- '
        )

class TestIndexFasta(unittest.TestCase):
    def setUp(self):
        self.seq_string = '\n'.join(['\n'.join(x) for x in aln])
        self.aln_stringio = StringIO(self.seq_string)

    def test_ensure_ordered_index(self):
        r = aln_pca.index_fasta(self.aln_stringio)
        for _aln, _seq in zip(aln,r):
            search_id = _aln[0].replace('>','')
            self.assertEqual(_aln[1], r[search_id])

class TestIdentityMatrix(unittest.TestCase):
    def setUp(self):
        self.seq_string = '\n'.join(['\n'.join(x) for x in aln])
        self.aln_stringio = StringIO(self.seq_string)
        self.index_fasta = aln_pca.index_fasta(self.aln_stringio)
    
    def test_buils_correct_matrix(self):
        r = aln_pca.identity_matrix(self.index_fasta)
        self.assertTrue(id_matrix.equals(r))
        self.assertTrue(id_matrix.columns.identical(r.columns))
        self.assertTrue(id_matrix.index.identical(r.index))
