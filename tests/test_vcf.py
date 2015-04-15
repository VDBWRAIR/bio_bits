import unittest
from bio_pieces import vcf_compare as vcfcat
from past.builtins import map
import re
import vcf


def fix_string(string):
    return re.sub(r'[\\\s\n\t]', '', string)


class TestVcfCompare(unittest.TestCase):
    def setUp(self):
       self.expected_head= '''AAQ  AC     ALT CB  CBD                                 CHROM   DP HPOLY  \
           0  NaN NaN  [None]  A  327  Den4/AY618992_1/Thailand/2001/Den4_1  327   NaN
       1  NaN NaN  [None]  G  325  Den4/AY618992_1/Thailand/2001/Den4_1  325   NaN
       2  NaN NaN  [None]  T  326  Den4/AY618992_1/Thailand/2001/Den4_1  326   NaN
       3  NaN NaN  [None]  T  326  Den4/AY618992_1/Thailand/2001/Den4_1  326   NaN
       4  NaN NaN  [None]  G  334  Den4/AY618992_1/Thailand/2001/Den4_1  334   NaN

          PAC  POS  PRC  QUAL  RAQ   RC REF
          0  NaN    1  100  None   39  327   A
          1  NaN    2  100  None   39  325   G
          2  NaN    3  100  None   39  326   T
          3  NaN    4  100  None   39  326   T
          4  NaN    5  100  None   39  334   G'''.strip()
       self.infile = open('tests/testinput/947.vcf')
       self.smallfile = 'tests/testinput/small.vcf'
       self.ambiguousfile = 'tests/testinput/ambigous.vcf'
       self.recs = list(vcf.Reader(open(self.smallfile)))
       self.ambiguous_recs = list(vcf.Reader(open(self.ambiguousfile)))
       self.smalldifferent = 'tests/testinput/smalldifferent.vcf'
       self.diff_recs = list(vcf.Reader(open(self.smalldifferent)))

    def get_positions(self, result):
        return [rec.POS for rec in result]

    def test_vcf_file_to_df(self):
        result_df = vcfcat.vcf_file_to_df(self.infile)
        raw_result_head = str(result_df.head()).strip()
        result_head, expected_head = map(fix_string, [raw_result_head, self.expected_head])
        self.assertEquals(expected_head, result_head)

    def test_filter_alt_equals(self):
        result = vcfcat.filter(self.recs, 'ALT', 'eq', ['C'])
        self.assertEquals(len(result) , 1)
        self.assertEquals(result[0].POS , 6)

    def test_filter_multi_alt_equals(self):
        result = vcfcat.filter(self.recs, 'ALT', 'eq', 'C, G')
        self.assertEquals(len(result) , 1)
        self.assertEquals(result[0].POS , 7)


    def test_filter_info_cbd(self):
        result = vcfcat.filter(self.recs, 'CBD', 'gt', 400)
        self.assertEquals(len(result) , 3)
        info = [(rec.POS, rec.REF) for rec in result]
        expected_info = [(8, 'A'), (9, 'G'), (10, 'T')]
        self.assertEquals(expected_info, info)


    def test_ambiguous_filter(self):
        result = vcfcat.ambiguous(self.ambiguous_recs)
        positions = self.get_positions(result)
        self.assertEquals([1, 8, 10], positions)
        cbs = [rec.INFO['CB'] for rec in result]
        self.assertEquals(['Y', 'R', 'Y'], cbs)


    def test_vcf_diff_str_no_threshold(self):
        result = vcfcat.diff(self.recs, self.diff_recs, 'REF', None)
        positions = self.get_positions(r[0] for r in result)
        self.assertEquals([1, 10], positions)

    def test_vcf_diff_int_with_threshold(self):
        result = vcfcat.diff(self.recs, self.diff_recs, 'DP', 100)
        positions = self.get_positions(r[0] for r in result)
        self.assertEquals([1, 3, 5], positions)

