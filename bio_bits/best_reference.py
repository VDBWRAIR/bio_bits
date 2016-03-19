from __future__ import division
from __future__ import print_function

from operator import itemgetter as get
from operator import attrgetter as attr
from itertools import groupby, starmap, chain, filterfalse
from typing import *
from functools import partial
#import sh
from plumbum.cmd import samtools
try:
    from functools import reduce
except:
    pass
"""
H3N2/CY074918/Managua/2010/NP_5
H3N2/CY074917/Managua/2010/NA_6
H3N2/CY074915/Managua/2010/HA_4
H3N2/CY074921/Managua/2010/PB1_2
H3N2/CY074920/Managua/2010/PA_3
H3N2/CY074919/Managua/2010/NS_8
H3N2/CY074916/Managua/2010/MP_7
H3N2/CY074922/Managua/2010/PB2_1
"""
segments = ["NP_5", "NA_6", "HA_4", "PB1_2", "PA_3", "NS_8", "MP_7", "PB2_1"]
#TODO:
# 1. for some reason `sh` wasn't working.
#    Need something that will do streams.
# 2. Need to aggregate the flu results to account for multiple segments
# 3. Encode this failure (and others):
#    [bam_idxstats] fail to load the index.

SeqID = str
Bam = str
OneBasedIndex = NamedTuple("OneBasedIndex", [("idx", int)])

IdxStatRow = NamedTuple("IdxStatRow",
                       [("refId", SeqID),
                        ("sequenceLength", int),
                        ("mappedReadCount", int),
                        ("unMappedReadCount", int)])

DepthRow = NamedTuple("DepthRow",
                    [("refId", SeqID),
                     ("position", OneBasedIndex),
                     ("depth", int)])

Results = NamedTuple("Results",
                   [("refId", SeqID),
                    ("sequenceLength", int),
                    ("coverageRatio", float),
                    ("avgDepth", float),
                    ('mappedRatio', float)])

DepthResults = NamedTuple("DepthResults",
                          [("refId", SeqID),
                           ("coverageRatio", float),
                           ("avgDepth", float)])

Options = NamedTuple("Options",
                     [("coverageWeight", Optional[float]),
                      ("mappedWeight", Optional[float]),
                      ("isFlu", bool)])

optionDefaults = {'coverageWeight' : 1.5, 'mappedWeight' : 1.0}

def makeDepthRow(values: List[str]) -> DepthRow:
   fields = ('refId', 'position', 'depth')
   d = dict(zip(fields, values))
   refId = SeqID(d['refId'])
   index = OneBasedIndex(int(d['position']))
   depth = int(d['depth'])
   return DepthRow(refId, index, depth)


def depth2row(dr: DepthRow) -> str:
    return '\t'.join([dr.refId, str(dr.position.idx), str(dr.depth)])

def idx2row(ir: IdxStatRow) -> str:
    return '\t'.join([ir.refId, str(ir.sequenceLength), str(ir.mappedReadCount), str(ir.unMappedReadCount)])

def makeIdxStatRow(values: List[str]) -> IdxStatRow:
   fields = ("refId", "sequenceLength", "mappedReadCount", "unMappedReadCount")
   d = dict(zip(fields, values))
   refId = SeqID(d['refId'])
   sequenceLength = int(d['sequenceLength'])
   mappedReadCount = int(d['mappedReadCount'])
   unMappedReadCount = int(d['unMappedReadCount'])
   return IdxStatRow(refId, sequenceLength, mappedReadCount, unMappedReadCount)

def getDepths(rawTSV: Iterable[str]) -> Iterable[DepthRow]:
    return map(makeDepthRow,  map(lambda x: x.split('\t'), filter(str.strip, rawTSV)))

def getIdxStats(rawTSV: Iterable[str]) -> Iterable[IdxStatRow]:
    # skip `*`, which is unmapped info, because it's not in Depth stats
    return filter(lambda x: x.refId != '*', map(makeIdxStatRow,  map(lambda x: x.split('\t'), filter(str.strip, rawTSV))))

def getDepthStats(rows: Iterable[DepthRow]) -> DepthResults:
    def accumulator(acc: Tuple[int,int,int], x: DepthRow) -> Tuple[int,int,int]:
        len, depthSum, uncovered = acc
        if x.depth == 0:
            return len + 1, depthSum, uncovered + 1
        else:
            return len+1, depthSum + x.depth, uncovered
    head = next(iter(rows))
    refId = head.refId
    rows = chain([head], rows)
    length, totalDepth, uncoveredCount = reduce(accumulator, rows, (0, 0, 0))
    coverageRatio =  (length - uncoveredCount) / length
    avgDepth = totalDepth / length
    return DepthResults(refId, coverageRatio, avgDepth)

def getMappedCount(row: IdxStatRow) -> int:
    return row.mappedReadCount

def makeResults(idxStat: IdxStatRow, depth: DepthResults) -> Results:
    assert depth.refId == idxStat.refId, "%s %s" % (depth, idxStat)
    return Results(idxStat.refId, idxStat.sequenceLength, depth.coverageRatio, depth.avgDepth, \
                   idxStat.mappedReadCount / (idxStat.mappedReadCount + idxStat.unMappedReadCount) )

Info = Union[IdxStatRow,DepthRow,DepthResults]

def allRefResults(idxStats: Iterable[IdxStatRow], depthStats: Iterable[DepthRow]) -> Iterable[Results]:
    # instead, get the results one at a time then match the results up
    getRef = lambda x: x.refId # type: Callable[[Info],SeqID]
    depthGroups = groupby(depthStats, key=getRef)
    depthInfo = map(getDepthStats, map(get(1), depthGroups))
    return map(makeResults, sorted(idxStats, key=getRef), sorted(depthInfo, key=getRef))

def samtoolsIdxStats(bam: Bam) -> Iterable[str]:
    return samtools['idxstats'][bam]().split('\n')
    #return sh.Command('samtools')('idxstats', bam, _iter=False)

def samtoolsDepth(bam: Bam) -> Iterable[str]:
    return samtools['depth'][bam]().split('\n')
     #return sh.Command('samtools')('depth', bam, _iter=False)

def sortByStats(xs: Iterable[Results], opts: Options) -> List[Results]:
    weighted = lambda x: (x.coverageRatio * opts.coverageWeight) + (x.mappedRatio * opts.mappedWeight) # type: Callable[[Results],float]
    return sorted(xs, key=weighted)

def find_sep(ref: str) -> str:
    seps = '/|:'
    for sep in seps:
         segmentsWithSeps = map(sep.__add__, segments)
         maybeSeps = list(filter(lambda seg: seg in ref, segmentsWithSeps))
         if maybeSeps:
             return maybeSeps[0][0]
    raise ValueError("seperator not found among %s in %s" % (seps, ref))

def combineResults(group: Tuple[SeqID, Iterable[Results]]) -> Results:
    ref, xs = group
    def combine(x: Results, y: Results) -> Results:
        return Results(ref,
                       x.sequenceLength + y.sequenceLength,
                       x.coverageRatio + y.coverageRatio,
                       x.avgDepth + y.avgDepth,
                       x.mappedRatio + y.mappedRatio)
    return reduce(combine, xs)

def avgFluSegments(resultsIter: Iterable[Results], opts: Options) -> List[Results]:
     results = list(resultsIter)
     sep = find_sep(results[0].refId)
     #justRef = lambda x: sep.join(filterfalse(segments.__contains__, x.refId.split(sep)))
     #TODO: this is a hack! the same references accross different segements don't necessarilyl have
     # the same ids even if you drop the segment from the string.
     justRef = lambda x: sep.join(list(filterfalse(segments.__contains__, x.refId.split(sep)))[2:])
     groups = groupby(sorted(results, key=justRef), key=justRef)
     combined = list(map(combineResults, groups))
     return combined

def getFinalStats(idxStats: Iterable[IdxStatRow], depthStats: Iterable[DepthRow], opts: Options) -> List[Results]:
    results = allRefResults(idxStats, depthStats)
    if opts.isFlu:
         results = avgFluSegments(results, opts)
    return sortByStats(results, opts)

def run(bam: Bam, opts: Options) -> None:
    idxStats  = getIdxStats(samtoolsIdxStats(bam))
    depthStats = getDepths(samtoolsDepth(bam))
    sortedStats = getFinalStats(idxStats, depthStats, opts)
    refList = '\n'.join(map(lambda x: x.refId, sortedStats))
    print(refList)

if __name__ == '__main__':
    from mypyextras.args import run_args
    import sys
    if sys.version[0] != '3':
        run.__annotations__ = dict(bam=Bam, opts=Options)
    run_args(run, optionDefaults)
