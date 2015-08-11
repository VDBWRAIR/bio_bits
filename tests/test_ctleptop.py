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
from bio_pieces.ctleptop import main, create_args
import subprocess
#from sys import version_info
# if version_info.major == 2:
    # import __builtin__ as builtins
# else:
    # import builtins

THIS = os.path.dirname(os.path.abspath(__file__))

class CtleptopFunctionTest(unittest.TestCase):

    """Base TestCase class, set up a CLI parser"""
    def setUp(self):
        self.infile = os.path.join(THIS, "Den4_MAAPS_TestData16.fasta")
        self.outfile = "out_file.fa"
        self.patch_fileinput = mock.patch('bio_pieces.ctleptop.access_mixed_aa')
        self.mock_fileinput = self.patch_fileinput.start()
        self.mock_fileinput_input = self.mock_fileinput.input
        self.addCleanup(self.mock_fileinput.stop)
        self.patch_sys_stdout = mock.patch('bio_pieces.ctleptop.sys.stdout')
        self.patch_sys_stderr = mock.patch('bio_pieces.ctleptop.sys.stderr')
        self.mock_stdout = self.patch_sys_stdout.start()
        self.mock_stderr = self.patch_sys_stderr.start()
        self.addCleanup(self.mock_stdout.stop)
        self.addCleanup(self.mock_stderr.stop)
        self.patcher_argparse = mock.patch('bio_pieces.ctleptop.create_args')
        self.mock_argparse = self.patcher_argparse.start()
        self.addCleanup(self.mock_argparse.stop)
        parser = ctleptop.create_args()
        self.parser = parser

    @mock.patch('bio_pieces.ctleptop.argparse')
    def test_args(self, mock_argparse):
        parser = ctleptop.create_args() #['-i' , self.infile, '-o', self.outfile])

        self.assertTrue(parser.i)
        self.assertTrue(parser.o)
        mock_argparse.ArgumentParser.return_value.add_argument.assert_call_once_with(
                ['-i' , self.infile, '-o', self.outfile]
        )

    #@mock.patch('bio_pieces.ctleptop.access_mixed_aa')
    #def test_access_mixed_aa(self, mock_mixed_aa):
        #ctleptop.access_mixed_aa(self.infile)
        #mock_mixed_aa.assert_called_with(self.infile)
        #my_mock = mock.MagicMock()
        #with mock.patch('__builtin__.open', my_mock):
            #manager = my_mock.return_value.__enter__.return_value
            #manager.read.return_value = ctleptop.access_mixed_aa(self.infile)
            #with open (self.infile) as h:
                #data = h.read()
            #my_mock.assert_called_once_with(self.infile)

    def test_main(self):
        cmd = ['ctleptop',
               '-i',
               self.infile,
               '-o',
               self.outfile,

        ]
        self.run_check_script(cmd)
    def run_check_script(self, cmd, stdin=None):
          if stdin is None:
              p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          else:
              p = subprocess.Popen(
                  cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                  stdin=stdin
              )
          sout, serr = p.communicate()
          self.assertEqual(b'Start processing and writing the output file to out_file.fa  please please wait ... \n',  sout)


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
