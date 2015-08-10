#!/usr/bin/env python
# encoding: utf-8
"""
test_ctleptop.py

Created by Dereje Jima on May 27, 2015
"""
import unittest
from bio_pieces import ctleptop
import mock
import tempfile
import os.path
from bio_pieces.ctleptop import main, parse_args
#from sys import version_info
# if version_info.major == 2:
    # import __builtin__ as builtins
# else:
    # import builtins

THISD = os.path.dirname(os.path.abspath(__file__))
print THISD




class CtleptopTestCase(unittest.TestCase):

    def test_with_empty_args(self):
        """User passes no args, should fail with SystemExit"""
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_sample_data(self):
        """ Te

        """
        args = self.parser.parse_args(['-i', self.example_fasta, '-o',
                                       self.outfile])
        result = ctl(args.tags)
        self.assertIsNotNone(result)


class CtleptopTest(unittest.TestCase):

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

    def test_isGap(self):
        result = ctleptop.isGap(["K/R", "I/T"], ["ARG", "NNN"])
        self.assertEquals(result, ["K/R", "GAPFOUND"])
