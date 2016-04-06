from bio_bits.best_reference import makeDepthRow, \
    makeIdxStatRow, depth2row, idx2row, allRefResults, \
    IdxStatRow, DepthRow, Results, OneBasedIndex, Options, \
    getFinalStats
from hypothesis import given, assume
import itertools
from hypothesis import strategies as st
from mypyextras.hypotypes import type_to_strat
import unittest

#TODO: test by negating the weight options; if they're negative, the order of results should be reversed.
#drAndIdx = st.lists(type_to_strat(DepthRow, {}), min_size=1).flatmap(
#    lambda xs:
#        map(lambda z: st.tuples(st.just(z), type_to_strat(IdxStatRow, {str: st.just(z.refId)}).filter(lambda x: x.mappedReadCount + x.unMappedReadCount > 0)), xs))
singleRefStats = type_to_strat(IdxStatRow, {}).filter(lambda x: x.mappedReadCount + x.unMappedReadCount > 0).flatmap(
    lambda x: st.tuples(st.just(x), st.lists(type_to_strat(DepthRow, {str: st.just(x.refId)}),min_size=1)))
#drAndIdx = st.lists(type_to_strat(DepthRow, {}).flatmap(
#    lambda x: st.tuples(st.just(z), type_to_strat(IdxStatRow, {str: st.just(z.refId)}).filter(lambda x: x.mappedReadCount + x.unMappedReadCount > 0)), xs))
#, min_size=1)
#@st.composite
#def drAndIdx(draw):
#    depths = st.lists(type_to_strat(DepthRow, {}), min_size=1).flatmap(
class TestBestReference(unittest.TestCase):

    def assertResultsEqual(self, x, y):
        self.assertEquals(x.refId, y.refId)
        self.assertEquals(x.sequenceLength, y.sequenceLength)
        self.assertAlmostEquals(x.coverageRatio, y.coverageRatio)
        self.assertAlmostEquals(x.mappedRatio, y.mappedRatio)
        self.assertAlmostEquals(x.avgDepth, y.avgDepth)

    @given(type_to_strat(DepthRow, {str : st.text().filter( lambda x: not ('\t' in x)) }))
    def test_depth_row_parse(self, row):
        self.assertEquals(row, makeDepthRow(depth2row(row).split('\t')))

    @given(type_to_strat(IdxStatRow, {str : st.text().filter( lambda x: not ('\t' in x)) }))
    def test_idx_parse(self, row):
        self.assertEquals(row, makeIdxStatRow(idx2row(row).split('\t')))


    def test_all_ref_results_example(self):
        depths = [DepthRow("foo", OneBasedIndex(12), 999), DepthRow("foo", OneBasedIndex(13), 12999)]
        idxs = [IdxStatRow("foo", 15, 44444, 90)]
        expected = [Results(refId='foo', sequenceLength=15,
                           coverageRatio=1.0, avgDepth=6999.0,
                           mappedRatio=0.9979790721695783)]
        actual = list(allRefResults(idxs, depths))
        self.assertEquals(len(expected), len(actual))
        for e, a in zip(expected, actual):
            self.assertResultsEqual(e, a)

    oneOrGreater = st.floats(min_value=1)

    def coverageRatio(depths):
        uncovered = list(filter(lambda x: x.depth == 0))
        return (len(depths) - len(uncovered)) / len(depths)

    mappedRatio = lambda x: x.mappedReadCount / (x.mappedReadCount + x.unMappedReadCount) 

    @given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def highest_mapped_count_ratio_is_first(self, rows, i):
        self.min_or_max_ratio_is_first_or_last(rows, i \
                ratio=MappedRatio,
                first=True)

    @given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def lowest_mapped_count_ratio_is_last(self, rows, i):
        self.min_or_max_ratio_is_first_or_last(rows, i \
                ratio=MappedRatio,
                first=False)

    @given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def highest_coverage_ratio_is_first(self, rows, i):
        self.min_or_max_ratio_is_first_or_last(rows, i \
                ratio=coverageRatio,
                first=True)

    @given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def lowest_coverage_ratio_is_last(self, rows, i):
        self.min_or_max_ratio_is_first_or_last(rows, i \
                ratio=coverageRatio,
                first=False)

    #@given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def min_or_max_ratio_is_first_or_last(self, rows, i, ratio, first):
        #if negative: i = -i 
        idxs, depths = zip(*rows)
        if str(ratio) == 'MappedRatio':
            opts = Options(coverageWeight=0,
                              mappedWeight=i, 
                              isFlu=False)
        else: 
            opts = Options(coverageWeight=i,
                              mappedWeight=0, 
                              isFlu=False)
        lastf = lambda x: list(x)[-1]
        firstf = lambda x: next(iter(x))
        if negative: 
            f, end = min, lastf
        else:
            f, end = max, firstf
        expectedEnd = f(idxs, key=ratio).refId
        res = getFinalStats(idxs, itertools.chain(*depths))
        acutalEnd = end(res).refId
        self.assertEquals(greatestMapped, expectedEnd) 


    #@given(st.lists(drAndIdx,min_size=1), posInt)
    @given(st.lists(singleRefStats,min_size=2), oneOrGreater)
    def test_negative_options_reverse_string_non_flu(self, rows, i):
        #TODO: This test fails because depths and idxs don't necessarily
        # have the same references
        #idx, depths = rows
        lmap = lambda f, xs: list(map(f, xs))
        idxs, depths = zip(*rows)
        ex_idx = lambda idx: (idx.sequenceLength, idx.mappedReadCount)#, idx.unMappedReadCount)
        ex_d = lambda d: d.depth
        idx_nums = lmap(ex_idx, idxs)
        d_nums = lmap(lambda ys: lmap(ex_d, ys), depths)
        # coverage ratio only cares about whether depth is non-zero
        # average depth is ignored
        # check that the results are not all equal
        def all_equal(xs):
            def _all_equal(xs, x):
                if not xs: return True
                else:
                    y = xs[0]
                    print(x, y)
                    return x == y and _all_equal(xs[1:], y)
            if not xs: raise ValueError("List was empty")
            else:
                return _all_equal(xs[1:], xs[0])
        assume(not all_equal(idx_nums))
        #assume(not(all_equal(idx_nums) and all_equal(d_nums)))

        depths = list(itertools.chain(*depths))
        #depths, idxs = zip(*rows)
        #depths, idxs = depths[0], [idxs[0]]
        posOpts = Options(coverageWeight=0,
                          mappedWeight=i, isFlu=False)
        posResult = getFinalStats(idxs, depths, posOpts)
        res_nums = lmap(lambda x:
                       (x.sequenceLength, x.coverageRatio,
                        x.mappedRatio), posResult)
        assume(not all_equal(res_nums))
        #assume(len(posResult) > 1)
        negOpts = Options(coverageWeight=0,
                          mappedWeight=i, isFlu=False)
        negResult = getFinalStats(idxs,depths, negOpts)

        for a,b in zip(posResult, reversed(negResult)):
            self.assertResultsEqual(a, b)

def test_extract_id_from_segment(self):
     ids = map(str.strip, """
     H3N2/CY074918/Managua/2010/NP_5
     H3N2/CY074917/Managua/2010/NA_6
     H3N2/CY074915/Managua/2010/HA_4
     H3N2/CY074921/Managua/2010/PB1_2
     H3N2/CY074920/Managua/2010/PA_3
     H3N2/CY074919/Managua/2010/NS_8
     H3N2/CY074916/Managua/2010/MP_7
     H3N2/CY074922/Managua/2010/PB2_1
     """.split('\n'))
     expected = "Managua/2010"
     for x in ids:
        self.assertEquals(expected, idFromSegmentFunc(x)(x))


#rows=[([DepthRow(refId='0', position=OneBasedIndex(idx=1), depth=1)],
#  IdxStatRow(refId='0', sequenceLength=1, mappedReadCount=1, unMappedReadCount=1))]
#, i=0)
