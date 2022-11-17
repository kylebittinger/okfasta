from pytest import approx

from okfastalib.msa import *


def test_msa_seqs():
    msa = MSA(["a", "b", "c"], ["ABC", "DEF", "GHI", "JKL"])
    assert list(msa.seqs) == [
        ("a", "ADGJ"), ("b", "BEHK"), ("c", "CFIL")
    ]

def test_msa_from_seqs():
    seqs = [("a", "ADGJ"), ("b", "BEHK"), ("c", "CFIL")]
    msa = MSA.from_seqs(seqs)
    assert msa.descs == ("a", "b", "c")
    assert msa.cols == ["ABC", "DEF", "GHI", "JKL"]

def test_msa_filter_by_index():
    msa = MSA(["a", "b", "c"], ["ABC", "DEF", "GHI", "JKL"])
    msa.filter_by_index((2, 4))
    assert msa.cols == ["DEF", "JKL"]
    
def test_msa_filter_by_index_remove():
    msa = MSA(["a", "b", "c"], ["ABC", "DEF", "GHI", "JKL"])
    msa.filter_by_index((2, 4), remove=True)
    assert msa.cols == ["ABC", "GHI"]

def test_msa_filter():
    msa = MSA(["a", "b", "c"], ["ABC", "DEF", "GHI", "JKL"])
    def has_no_Gs(xs):
        return "G" not in xs
    msa.filter(has_no_Gs)
    assert msa.cols == ["ABC", "DEF", "JKL"]

def test_msa_map():
    msa = MSA(["a", "b", "c"], ["ABC", "DEF", "GHI", "JKL"])
    def make_lowercase(xs):
        return xs.lower()
    assert list(msa.map(make_lowercase)) == ["abc", "def", "ghi", "jkl"]

def test_strcat():
    assert strcat(['ab', 'cde']) == 'abcde'

def test_shannon():
    assert shannon([31, 0, 0, 0]) == 0
    assert shannon([5, 5, 0, 0]) == approx(0.6931472)

def test_column_stats():
    m = MSA('abcdefghij', ['aaaabbbb--'])
    s = next(m.column_stats())
    assert s == {
        "column_position": 1,
        "number_of_values": 8,
        "gaps_proportion": 0.2,
        "entropy": approx(0.6931472),
        "consensus_value": "a",
        "consensus_proportion": 0.5
    }

def test_column_stats_allgaps():
    m = MSA('abcdefghij', ['----------'])
    s = next(m.column_stats())
    assert s == {
        "column_position": 1,
        "number_of_values": 0,
        "gaps_proportion": 1.0,
        "entropy": 0.0,
        "consensus_value": "-",
        "consensus_proportion": 1.0
    }

def test_column_stats_nogaps():
    m = MSA('abcdefghij', ['ccggttaacc'])
    s = next(m.column_stats())
    assert s == {
        "column_position": 1,
        "number_of_values": 10,
        "gaps_proportion": 0.0,
        "entropy": approx(1.332179),
        "consensus_value": "c",
        "consensus_proportion": 0.4
    }
