import unittest
from bio_pieces import vcf_compare as mod
from past.builtins import map

def fix_string(string):
    import re
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
       self.infile = 'tests/testinput/947.vcf'


    def test_vcf_file_to_df(self):
        result_df = mod.vcf_file_to_df(self.infile)
        raw_result_head = str(result_df.head()).strip()
        result_head, expected_head = map(fix_string, [raw_result_head, self.expected_head])
        self.assertEquals(expected_head, result_head)







