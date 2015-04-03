from os.path import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest

import mock

from bio_pieces import amos

RED =   '{{RED\n' \
        'iid:{iid}\n' \
        'eid:{iid}\n' \
        'seq:\n' \
        '{seq}\n' \
        '.\n' \
        'qlt:\n' \
        '{qlt}\n' \
        '.\n' \
        '}}\n'
CTG =   '{{CTG\n' \
        'iid:{iid}\n' \
        'eid:{eid}\n' \
        'com:\n' \
        '{com}\n' \
        '.\n' \
        'seq:\n' \
        '{seq}\n' \
        '.\n' \
        'qlt:\n' \
        '{qlt}\n' \
        '.\n' \
        '{tlelist}' \
        '}}\n'
TLE =   '{{TLE\n' \
        'src:{src}\n' \
        'off:{off}\n' \
        'clr:{clr}\n' \
        '}}\n'
AMOS =  '{redlist}{ctglist}'

def make_red(iid, eid, seq, qlt):
    return RED.format(
        iid=iid, eid=eid, seq=seq, qlt=qlt
    )

def make_tle(src, off, clr):
    return TLE.format(
        src=src, off=off, clr=clr
    )

def make_ctg(iid, eid, com, seq, qlt, tlelist):
    return CTG.format(
        iid=iid, eid=eid, com=com, seq=seq, qlt=qlt,
        tlelist=''.join(tlelist)
    )

class TestSplitLine(unittest.TestCase):
    def test_splits_int(self):
        self.assertEqual(('a',1), amos.splitline('a:1', ':', int))
    
    def test_splits_float(self):
        self.assertEqual(('a',1.0), amos.splitline('a:1', ':', float))

    def test_raises_exception_not_enough_pieces(self):
        self.assertRaises(ValueError, amos.splitline, 'a', ':', int)

    def test_raises_exception_too_many_pieces(self):
        self.assertRaises(ValueError, amos.splitline, 'a:1:1', ':', int)

class TestAmos(unittest.TestCase):
    def setUp(self):
        self.redlist = []
        self.tlelist = []
        self.ctglist = []
        for i in range(10):
            redstr = make_red(i, i, 'ATGC', 'DDDD')
            tlestr = make_tle(i, 0, '0,4')
            self.redlist.append(redstr)
            self.tlelist.append(tlestr)

        for i in range(5):
            ctgstr = make_ctg(
                i, 'foo-{0}'.format(i), 'foocom', 'ATGC', 'DDDD', self.tlelist[i:i+2]
            )
            self.ctglist.append(ctgstr)

        self.amosstr = ''.join(self.redlist) + ''.join(self.ctglist)

    def test_parses_correct_amos_contents(self):
        mock_fh = mock.MagicMock()
        mock_fh.__iter__.return_value = self.amosstr.splitlines(True)
        _amos = amos.AMOS(mock_fh)
        self.assertEqual(10, len(_amos.reds))
        self.assertEqual(5, len(_amos.ctgs))
        self.assertIsInstance(_amos.reds, dict)
        self.assertIsInstance(_amos.reds[0], amos.RED)
        self.assertIsInstance(_amos.ctgs, dict)
        self.assertIsInstance(_amos.ctgs[0], amos.CTG)

class TestCtg(unittest.TestCase):
    def setUp(self):
        self.tlelist = []
        self.reddict = {}
        for i in range(10):
            red_str = make_red(i,i,'ATGC','DDDD')
            red = amos.RED.parse(red_str)
            self.reddict[red.iid] = red
            self.tlelist.append(
                make_tle(i,'0','0,4')
            )
        self.ctg_str = make_ctg(1,'foobar','foo','ATGC','DDDD', self.tlelist)

    def test_parses_valid_ctg_block(self):
        ctg = amos.CTG.parse(self.ctg_str)
        self.assertEqual(1, ctg.iid)
        self.assertEqual('foobar', ctg.eid)
        self.assertEqual('foo', ctg.com)
        self.assertEqual('ATGC', ctg.seq)
        self.assertEqual('DDDD', ctg.qlt)
        self.assertEqual(10, len(ctg.tlelist))
        firsttle = ctg.tlelist[0]
        lasttle = ctg.tlelist[-1]
        self.assertEqual(0, firsttle.src)
        self.assertEqual(9, lasttle.src)

    def test_parses_ctg_with_no_tle(self):
        ctgstr = make_ctg(1,'foobar','foo','ATGC','DDDD', [])
        ctg = amos.CTG.parse(ctgstr)
        self.assertEqual(0, len(ctg.tlelist))

class TestTle(unittest.TestCase):
    def setUp(self):
        self.tle_str = make_tle(1, 0, '0,4')

    def test_parses_valid_tle_block(self):
        tle = amos.TLE.parse(self.tle_str)
        self.assertEqual(1, tle.src)
        self.assertEqual(0, tle.off)
        self.assertEqual('0,4', tle.clr)

    def test_raises_exception_with_bad_tle_block(self):
        self.assertRaises(ValueError, amos.TLE.parse, '')

class TestRed(unittest.TestCase):
    def setUp(self):
        self.red_str = make_red(1,1,'ATGC','DDDD')
        self.seqrec = mock.Mock(
            id='read1',
            seq='ATGC',
            _per_letter_annotations={'phred_quality':'IIII'}
        )

    def test_returns_red_instance(self):
        self.assertIsInstance(amos.RED.parse(self.red_str), amos.RED)
        
    def test_parses_valid_red_block(self):
        red = amos.RED.parse(self.red_str)
        self.assertEqual(1, red.iid)
        self.assertEqual(1, red.eid)
        self.assertEqual('ATGC', red.seq)
        self.assertEqual('DDDD', red.qlt)

    def test_raises_exception_with_bad_red_block(self):
        self.assertRaises(ValueError, amos.RED.parse, '')

    def test_sets_from_seqrecord(self):
        red = amos.RED.parse(self.red_str)
        red.set_from_seqrec(self.seqrec)
        self.assertEqual('ATGC', red.seq)
        self.assertEqual('IIII', red.qlt)
        self.assertEqual('read1', red.eid)

    def test_raises_exception_with_mismatch_seqrec_seq(self):
        red = amos.RED.parse(self.red_str)
        self.seqrec.seq = 'AAAA'
        self.assertRaises(
            amos.RED.MismatchedSequenceError,
            red.set_from_seqrec, self.seqrec
        )

    def test_raises_exception_with_missing_phred_quality(self):
        red = amos.RED.parse(self.red_str)
        del self.seqrec._per_letter_annotations['phred_quality']
        self.assertRaises(
            ValueError,
            red.set_from_seqrec, self.seqrec
        )