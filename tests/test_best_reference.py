from bio_bits.best_reference import makeDepthRow, \
    makeIdxStatRow, depth2row, idx2row, allRefResults, \
    IdxStatRow, DepthRow, Results, OneBasedIndex, Options, \
    getFinalStats
from hypothesis import given
from hypothesis import strategies as st
from mypyextras.hypotypes import type_to_strat
import unittest

#TODO: test by negating the weight options; if they're negative, the order of results should be reversed.
drAndIdx = st.lists(type_to_strat(DepthRow, {}), min_size=1).flatmap(
    lambda x: st.tuples(st.just(x), type_to_strat(IdxStatRow, {})))
posInt = st.integers(min_value=0)
class TestBestReference(unittest.TestCase):

    @given(type_to_strat(DepthRow, {str : lambda x: not ('\t' in x) }))
    def test_depth_row_parse(self, row):
        self.assertEquals(row, makeDepthRow(depth2row(row).split('\t')))

    @given(type_to_strat(IdxStatRow, {str : lambda x: not ('\t' in x) }))
    def test_idx_parse(self, row):
        self.assertEquals(row, makeIdxStatRow(idx2row(row).split('\t')))

    def assertResultsEqual(self, x, y):
        self.assertEquals(x.refId, y.refId)
        self.assertEquals(x.sequenceLength, y.sequenceLength)
        self.assertAlmostEquals(x.coverageRatio, y.coverageRatio)
        self.assertAlmostEquals(x.mappedRatio, y.mappedRatio)
        self.assertAlmostEquals(x.avgDepth, y.avgDepth)

    def test_all_ref_results_example(self):
        depths = [DepthRow("foo", OneBasedIndex(12), 999), DepthRow("foo", OneBasedIndex(13), 12999)]
        idxs = [IdxStatRow("foo", 15, 44444, 90)]
        expected = Results(refId='foo', sequenceLength=15,
                           coverageRatio=1.0, avgDepth=6999.0,
                           mappedRatio=0.9979790721695783)
        actual = next(iter(allRefResults(idxs, depths)))
        self.assertResultsEqual(expected, actual)

#    @given(st.lists(drAndIdx,min_size=1), posInt)
#    def test_negative_options_reverse_string_non_flu(self, rows, i):
#        #TODO: This test fails because depths and idxs don't necessarily
#        # have the same references
#        depths, idxs = zip(*rows)
#        depths, idxs = depths[0], [idxs[0]]
#        posOpts = Options(i, 1, False)
#        posResult = getFinalStats(idxs, depths, posOpts)
#        negOpts = Options(-i, 1, False)
#        negResult = getFinalStats(idxs, depths, negOpts)
#        print(negResult, posResult)
#        for a,b in zip(posResult, reversed(negResult)):
#            self.assertResultsEqual(a, b)



#rows=[([DepthRow(refId='0', position=OneBasedIndex(idx=1), depth=1)],
#  IdxStatRow(refId='0', sequenceLength=1, mappedReadCount=1, unMappedReadCount=1))]
#, i=0)
