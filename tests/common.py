import os
import sys
from os.path import *
import shutil
import tempfile
import subprocess

from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr
from mock import Mock, patch
import tempdir

# tests dir
TESTDIR = dirname(abspath(__file__))
EXAMPLES = join(TESTDIR, 'example_files')
 
class BaseTester(object):
    def setUp(self):
        ''' Auto enter temporary directory '''
        self.tdir = tempdir.TempDir()
        os.chdir(self.tdir.name)
     
    def tearDown(self):
        ''' Clean up the tempdir after test '''
        self.tdir.dissolve()
     
    def _C( self, *args, **kwargs):
        '''
        Set modulepath as instance variable to be the name of the module
        Set functionname as instance variable to automagically run that function
        with self._C
        '''
        return self._run_func( self.modulepath, self.functionname, *args, **kwargs )

    def _run_func( self, modulepath, functionname, *args, **kwargs ):
        ''' Run a function from modulepath '''
        m = __import__( modulepath, fromlist=[self.functionname] )
        return getattr(m,functionname)( *args, **kwargs )

    def make_seqrec(self, seq, quals, id='id'):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        seq = Seq( seq, generic_dna )
        rec = SeqRecord(
            seq,
            id=id,
            description='description',
            name='name'
        )
        rec._per_letter_annotations['phred_quality'] = quals
        return rec

    def random_seqs(self, numseqs=100):
        import random
        dna = 'ATGC'
        seqs = []
        maxlen = 0.0
        maxqual = 0.0
        for i in range( 1, numseqs ):
            randbase = random.choice( dna )
            randlen = random.randint( 0, 1000 )
            randseq = ''.join( [random.choice(dna) for i in range(0,randlen)] )
            randqual = [random.randint(0,60) for i in range(0, randlen)]
            aqual = 0
            if randlen != 0:
                aqual = round( sum(randqual) * 1.0 / randlen )
                maxqual = max( aqual, maxqual )
            maxlen = max( maxlen, randlen )
            seqs.append( make_seqrec( randseq, randqual ) )
        return (seqs, maxlen, maxqual)

    def rand_seqrec(self, seqlen=10, cal=0, car=0, cql=0, cqr=0):
        '''
            Makes a Bio.SeqRecord.SeqRecord such that the .seq returns
            a Bio.Seq.Seq with a random sequence that has length seqlen if seqlen is int otherwise makes Seq with seqlen
            Then the SeqRecord has ._per_letter_annotations['phred_quality'] with random qualities
            for all the bases
            SeqRecord.annotations clip_adapter_left, right, clip_qual_left, right are set to cal, car, cql and cqr
        '''
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        import random
        seq = None
        if isinstance( seqlen, int ):
            # Random Sequence
            seq = Seq( rand_seq(seqlen), generic_dna )
            # Random qualities
            qual = [random.randint(1,40) for i in range(seqlen)]
        else:
            seq = Seq( seqlen, generic_dna )
            # Random qualities
            qual = [random.randint(1,40) for i in range(len(seqlen))]
        # Random id. Hopefully random enough so no duplicate ids
        id = 'seq_{}'.format(random.randint(1,999999999999))
        record = SeqRecord(
            seq,
            id=id,
            description='Random sequence',
            name=id
        )
        record._per_letter_annotations['phred_quality'] = qual
        record.annotations['clip_adapter_left'] = cal
        record.annotations['clip_adapter_right'] = car
        record.annotations['clip_qual_left'] = cql
        record.annotations['clip_qual_right'] = cqr
        return record

    def rand_quals(self, qlen=10, maxq=40):
        '''return random quals'''
        import random
        return [random.randint(0,maxq) for i in range(qlen)]

    def rand_seq(self, seqlen=10):
        ''' return random sequence length seqlen '''
        import random
        dna = ('A','C','G','T')
        return ''.join( [dna[random.randint(0,3)] for i in range(seqlen)] )
