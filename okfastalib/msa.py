from dataclasses import dataclass
import collections
import math
import itertools


def enumerate1(xs):
    for n, x in enumerate(xs):
        yield (n + 1), x

def range1(stop):
    for n in range(stop):
        yield (n + 1)

class MSA:
    def __init__(self, descs, cols):
        self.descs = descs
        self.cols = cols

        # Double check that length of each column is equal to number
        # of sequence descriptions
        len_descs = len(descs)
        for col in self.cols:
            assert len(col) == len_descs

    def filter(self, fcn):
        self.cols = [col for col in self.cols if fcn(col)]
        return self

    def filter_by_index(self, idxs, remove=False):
        idxs = set(idxs)
        if remove:
            idxs = [idx for idx in range1(len(self.cols)) if idx not in idxs]
        self.cols = [col for idx, col in enumerate1(self.cols) if idx in idxs]
        return self

    def map(self, fcn):
        return map(fcn, self.cols)

    @property
    def seqs(self):
        seqvals = map(strcat, zip(*self.cols))
        yield from zip(self.descs, seqvals)

    @classmethod
    def from_seqs(cls, seqs):
        descs, seqvals = list(zip(*seqs))
        cols = [
            strcat(col_chars) for col_chars
            in itertools.zip_longest(*seqvals, fillvalue="-")
        ]
        return cls(descs, cols)

@dataclass
class ColumnStatsResult:
    nvals: int
    prop_gaps: float
    entropy: float
    consensus_val: str
    consensus_prop: float

def column_stats(col):
    len_col = len(col)
    ngaps = col.count("-")
    nvals = len_col - ngaps
    if nvals == 0:
        return ColumnStatsResult(
            nvals = 0,
            prop_gaps = 1.0,
            entropy = 0.0,
            consensus_val = "-",
            consensus_prop = 1.0,
        )
    assert nvals > 0
    ctr = collections.Counter(col)
    del ctr["-"]
    cts = ctr.values()
    consensus_val, consensus_cts = ctr.most_common(1)[0]
    return ColumnStatsResult(
        nvals = nvals,
        prop_gaps = ngaps / len_col,
        entropy = shannon(cts),
        consensus_val = consensus_val,
        consensus_prop = consensus_cts / nvals,
    )

def shannon(cts):
    total = sum(cts)
    props = [c / total for c in cts]
    shannon_terms = [p * math.log(p) if p > 0 else 0 for p in props]
    return -sum(shannon_terms)

def strcat(xs):
    return ''.join(xs)
