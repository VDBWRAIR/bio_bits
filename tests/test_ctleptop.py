#!/usr/bin/env python
# encoding: utf-8
"""
test_ctleptop.py

Created by Dereje Jima on May 27, 2015
"""
import unittest
import os
from bio_pieces import ctleptop as ctl
from mock import patch
#from sys import version_info
# if version_info.major == 2:
    # import __builtin__ as builtins
# else:
    # import builtins

THISD = os.path.dirname(os.path.abspath(__file__))


# class CommandLineTestCase(unittest.TestCase):

#     """Base TestCase class, setup a CLI parser """

#     def setUpClass(self):
#         parser = ctl.parse_args()
#         ctl.parser = parser
#         self.example_fasta = os.path.join(THISD, '721TestConsensus.fasta')
#         self.outfile = os.path.join(THISD, "outfile.fa")

#     @patch(ctl.main)
#     def test_main(self, mockopt):
#         mockopt.return_value = {'-i': self.example_fasta, '-o': self.outfile}
#         lambda fn, fn2: (open(fn).read(), open(fn2).read())
#         ctl.main()
#         map(self.assertTrue, map(os.path.exists, self.outfile))


# class CtleptopTestCase(CommandLineTestCase):

#     def test_with_empty_args(self):
#         """User passes no args, should fail with SystemExit"""
#         with self.assertRaises(SystemExit):
#             self.parser.parse_args([])

#     def test_sample_data(self):
#         """ Te

#         """
#         args = self.parser.parse_args(['-i', self.example_fasta, '-o',
#                                        self.outfile])
#         result = ctl(args.tags)
#         self.assertIsNotNone(result)


class CtleptopTest(unittest.TestCase):

    """Docstring for CtleptopTest. """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getaalist(self):
        result = ctl.getaalist(['AAA', 'ACT'])
        self.assertEquals(result, ['K', 'T'])

    def test_nearbyPermutations(self):
        result = ctl.nearbyPermutations("AAR")
        self.assertEquals(result, set(['AAG', 'AAA']))
        result2 = ctl.nearbyPermutations("ARR")
        self.assertEquals(result2, set(['AGG', 'AAG', 'AAA', 'AGA']))
        result3 = ctl.nearbyPermutations("AAA")
        self.assertEquals(result3, set(['AAA']))

    def test_permutations(self):
        """Test listReplace function.
        :returns: bool

        """
        result = ctl.permutations(set(['CA']), ['A', 'T'])
        self.assertEquals(result, set(['ACA', 'TCA']))

    def test_list_overlap(self):
        """Check if two list overlap
        :returns: bool

        """
        result = ctl.list_overlap(
            "GAY", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result)
        result2 = ctl.list_overlap(
            "WTA", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result2)
        result3 = ctl.list_overlap(
            "TTT", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(False,  result3)

    def test_isGap(self):
        result = ctl.isGap(["K/R", "I/T"], ["ARG", "NNN"])
        self.assertEquals(result, ["K/R", "GAPFOUND"])
