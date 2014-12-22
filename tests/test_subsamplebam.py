from common import *

from bio_pieces import subsamplebam

class Base(BaseTester):
    # mpileup template without mapping quals(no -s option to samtools mpileup)
    PILE_TEMPLATE = '{ref}\t{pos}\t{refbase}\t{depth}\t{bases}\t{quals}\n'
    # mpileup template with mapping quals(with -s option to samtools mpileup)
    PILE_TEMPLATE_MQUALS = PILE_TEMPLATE[:-1] + '\t{mquals}\n'

    def _mpileup_factory(self, mpileups):
        '''
        Turn a bunch of mpileup strings into a mocked mpileup call
        '''
        def _mock_mpileup(self, *args):
            '''
            Mock out the mpileup call
            '''
            return mpileups

    def _format_mpileup(self, *args, **kwargs):
        '''
        Make mpileup from PILE_TEMPLATE
        '''
        if kwargs.get('ref',False):
            if kwargs.get('mquals',False):
                return self.PILE_TEMPLATE_MQUALS.format(**kwargs)
            else:
                return self.PILE_TEMPLATE.format(**kwargs)
        else:
            if len(args) == 6:
                return self.PILE_TEMPLATE_MQUALS.format(**kwargs)
            else:
                return self.PILE_TEMPLATE.format(**kwargs)

@patch('bio_pieces.sbusamplebam.Popen')
class TestMpileup(Base):
    def setUp(self):
        super(TestMpileup,self).setUp()

        self.bam = 'test.bam'
        self.refname = 'Ref1'
        self.reffile = 'test.fasta'
        self.refseq = 'ATGC'*2 + 'AT'
        self.restart = 1
        self.refend = len(self.refseq)
        self.refrecord = Mock(seq=self.refseq, id=self.refname)
        # -q
        self.map_q = 0
        # -Q
        self.base_q = 0
        # -d
        self.max_depth = 1000
        # -r
        self.regionstring = '{0}:{1}-{2}'.format(self.refname, self.refstart, self.refend)
        self.mpileup_args = [
            '-q', self.map_q, '-Q', self.base_q, '-d', self.max_depth, '-s', '-r', self.regionstring, self.bam
        ]

    def test_uses_regionstring(self, mpopen):

    def test_uses_dash_s_mapping_quality(self, mpopen):
        pass

    def test_limits_based_on_params(self, mpopen):
        pass

    def test_using_f_sets_refbase(self, mpopen):

    def test_no_f_or_s(self, mpopen):
        # No f means all refbase = N
        # No s means no mapping quals
        pileups = []
        for i in range(1,11):
            self._format_pileup(self.refname, i, 'N', 10, 'A'*10, 'I'*10)
        mpopen.return_value.stdout = iter(pileups)
