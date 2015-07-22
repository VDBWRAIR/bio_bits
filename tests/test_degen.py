import unittest
import mock
from bio_pieces.degen import Gene, get_gene_pos_seq, get_result_table
from itertools import starmap


class DegenTest(unittest.TestCase):
    def setUp(self):
        self.seq = 'G'*83 + 'GGRAAY' + 420*'A' + 'RAGT'+'C'*8000 + 'YYY' +  'A' * 2000
        self.positions = [85, 88, 509, 8513, 8514, 8515]
        self.genestrs = ["anchored capsid protein"] * 2 + ["membrane glycoprotein precursor"] + ["nonstructural protein NS5"]*3
        self.nts = [self.seq[i] for i in self.positions]
        self.genbank_id = 'KJ189367.1'
        self.expected_str =''
        for p, g, nt in zip(self.positions, self.genestrs, self.nts):
            self.expected_str += '\t'.join( map(str, [g, p, nt]) ) + '\n'
        self.expected_str = self.expected_str[:-1]
        self.geneinfo =    [ ["anchored capsid protein",             84, 426],
        ["membrane glycoprotein precursor", 426,923,     ],
        ["envelope protein",                924,2408,    ],
        ["nonstructural protein NS1",       2409,3464,   ],
        ["nonstructural protein NS2A",      3465,4118,   ],
        ["nonstructural protein NS2B",      4119,4508,   ],
        ["nonstructural protein NS3",       4509,6365,   ],
        ["nonstructural protein NS4A",      6366,6746,   ],
        ["2K peptide",                      6747,6815,   ],
        ["nonstructural protein NS4B",      6816,7562,  ],
        ["nonstructural protein NS5",        7563,10259,]]
        self.genes = starmap(Gene, self.geneinfo)


    @mock.patch('Bio.SeqIO.parse')
    def test_functional(self, mparse):
        '''Note: requires internet access to genbank.'''
        mparse.return_value = [self.seq]
        actual = next(get_result_table('_', ref_id=self.genbank_id))
        self.assertMultiLineEqual(actual, self.expected_str)

    def test_get_gene_pos_seq(self):
        expected = sorted(zip(self.genestrs, self.positions, self.nts))
        actual = sorted(get_gene_pos_seq(self.genes, self.seq))
        self.assertListEqual(actual, expected)


