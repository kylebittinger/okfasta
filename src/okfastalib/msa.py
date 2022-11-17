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

    def filter_by_index(self, idxs, remove=False):
        idxs = set(idxs)
        if remove:
            idxs = [idx for idx in range1(len(self.cols)) if idx not in idxs]
        self.cols = [col for idx, col in enumerate1(self.cols) if idx in idxs]
        return self

    column_stats_header = [
            "column_position", "number_of_values", "gaps_proportion",
            "entropy", "consensus_value", "consensus_proportion",
        ]

    column_stats_fmt = "{0}\t{1}\t{2:1.2f}\t{3:1.4f}\t{4}\t{5:1.2f}"

    def column_stats(self):
        for col_position, col in enumerate1(self.cols):
            len_col = len(col)
            ngaps = col.count("-")
            nvals = len_col - ngaps
            if nvals == 0:
                yield {
                    "column_position": col_position,
                    "number_of_values": 0,
                    "gaps_proportion": 1.0,
                    "entropy": 0.0,
                    "consensus_value": "-",
                    "consensus_proportion": 1.0
                }
                continue
            ctr = collections.Counter(col)
            del ctr["-"]
            cts = ctr.values()
            consensus_val, consensus_cts = ctr.most_common(1)[0]
            yield {
                "column_position": col_position,
                "number_of_values": nvals,
                "gaps_proportion": ngaps / len_col,
                "entropy": shannon(cts),
                "consensus_value": consensus_val,
                "consensus_proportion": consensus_cts / nvals,
            }

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


def shannon(cts):
    cts = [c for c in cts if c > 0]
    # If we use the formula when h=0, python will return -0.0
    if len(cts) == 1:
        return 0.0
    total = sum(cts)
    props = [c / total for c in cts]
    h = -sum(p * math.log(p) for p in props)
    return h


def strcat(xs):
    return ''.join(xs)
