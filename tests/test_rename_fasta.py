try:
    import unittest2 as unittest
except ImportError:
    import unittest

import shutil
import subprocess
import os
import StringIO
from os.path import *

import mock

from bio_pieces import rename_fasta

THIS = dirname(__file__)
TEST_FASTA = join(THIS,'example.fasta')
TEST_CSV = join(THIS,'renamelist.csv')

class TestGetCSVMap(unittest.TestCase):
    def test_gets_mapping(self):
        mapping = rename_fasta.get_csv_map(TEST_CSV)
        self.assertEqual(mapping['001'], 'sample1')
        self.assertEqual(mapping['002'], 'sample2')
        #self.assertEqual(mapping['003'], 'sample3')
        self.assertNotIn('#From', mapping)
        self.assertNotIn('003', mapping)

class TestFunctional(unittest.TestCase):
    def setUp(self):
        self.test_fasta = '/tmp/test.fasta'
        shutil.copyfile(TEST_FASTA, self.test_fasta)
        self.addCleanup(os.unlink, self.test_fasta)
        self.mapping = rename_fasta.get_csv_map(TEST_CSV)
        self.rev_mapping = {}
        for key,value in self.mapping.items():
            self.rev_mapping[value] = key
        self.rev_mapping['003'] = '003'

        self.patcher_argparse = mock.patch('bio_pieces.rename_fasta.argparse')
        self.mock_argparse = self.patcher_argparse.start()
        self.addCleanup(self.mock_argparse.stop)

    def test_renames_inplace(self):
        cmd = [
            'rename_fasta',
            TEST_CSV,
            self.test_fasta,
            '--inplace'
        ]

        self._run_check_script(cmd, self.rev_mapping)

    def test_renames_stdin_stdout(self):
        cmd = [
            'rename_fasta',
            TEST_CSV,
            '-',
        ]

        sin = open(self.test_fasta)

        self._run_check_script(cmd, self.rev_mapping, sin)

    def _run_check_script(self, cmd, rev_mapping, stdin=None):
        if stdin is None:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                stdin=stdin
            )
        sout, serr = p.communicate()

        if '--inplace' in cmd:
            _iter = open(self.test_fasta)
        else:
            _iter = sout.splitlines()

        for line in _iter:
            if line.startswith('>'):
                p = line.split()
                id = p[0][1:]
                self.assertIn(id, rev_mapping)

        self.assertEqual('003 is not in provided mapping\n', serr)
