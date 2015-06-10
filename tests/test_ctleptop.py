#!/usr/bin/env python
# encoding: utf-8
"""
test_ctleptop.py

Created by Dereje Jima on May 27, 2015
"""
import unittest
import os
import sys
from bio_pieces import ctleptop as ctl

THISD = os.path.dirname(os.path.abspath(__file__))
class CtleptopTest(unittest.TestCase):

    """Docstring for CtleptopTest. """

    def setUp(self):
        self.example_fasta=os.path.join(THISD, '721TestConsensus.fasta')

    def tearDown(self):
        pass

    def test_readFasta(self):
        pass
    def test_listReplace(self):
        """Test listReplace function.
        :returns: bool

        """
        result = ctl.listReplace("GAY", "Y", ["C", "T"])
        self.assertEquals(result, ['D', 'D'])
        result2 = ctl.listReplace("ARA", "R", ["A", "G"])
        self.assertEquals(result2, ['K', 'R'])
        result3 = ctl.listReplace("WTA", "W", ["A", "T"])
        self.assertEquals(result3, ['I', 'L'])

    def test_list_overlap(self):
        """Check if two list overlap
        :returns: bool

        """
        result = ctl.list_overlap("GAY", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result)
        result2 = ctl.list_overlap("WTA", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result2)
        result3 = ctl.list_overlap("TTT", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(False,  result3)



