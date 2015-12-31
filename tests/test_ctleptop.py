import mock
import os

from bio_bits.compat import unittest, StringIO

from bio_bits import degen
from bio_bits import ctleptop

class TestCtleptopTest(unittest.TestCase):
    """Docstring for CtleptopTest. """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getaalist(self):
        result = ctleptop.getaalist(['AAA', 'ACT'])
        self.assertEquals(result, ['K', 'T'])

    def test_nearbyPermutations(self):
        result = ctleptop.nearbyPermutations("AAR")
        self.assertEquals(result, set(['AAG', 'AAA']))
        result2 = ctleptop.nearbyPermutations("ARR")
        self.assertEquals(result2, set(['AGG', 'AAG', 'AAA', 'AGA']))
        result3 = ctleptop.nearbyPermutations("AAA")
        self.assertEquals(result3, set(['AAA']))

    def test_permutations(self):
        """Test listReplace function.
        :returns: bool

        """
        result = ctleptop.permutations(set(['CA']), ['A', 'T'])
        self.assertEquals(result, set(['ACA', 'TCA']))

    def test_list_overlap(self):
        """Check if two list overlap
        :returns: bool

        """
        result = ctleptop.list_overlap(
            "GAY", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result)
        result2 = ctleptop.list_overlap(
            "WTA", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result2)
        result3 = ctleptop.list_overlap(
            "TTT", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(False,  result3)

class TestModEntry(unittest.TestCase):
    def setUp(self):
        self.entry = ('SEQUENCE', 9, 3, 'AAA', 'R', 'GENE')

    def test_returns_unmod_entry(self):
        r = ctleptop.mod_entry(self.entry, degen.Gene('CDS', 1, 10))
        self.assertEqual(self.entry, r)

    def test_returns_noncoding_entry(self):
        r = ctleptop.mod_entry(self.entry, degen.Gene('CDS', 1, 1))
        self.assertEqual(r[4], 'NON-CODING')

    def test_returns_gapfound_entry(self):
        self.entry = list(self.entry)
        self.entry[3] = 'NNN'
        self.entry = tuple(self.entry)
        r = ctleptop.mod_entry(self.entry, degen.Gene('CDS', 1, 10))
        self.assertEqual(r[4], 'GAPFOUND')
