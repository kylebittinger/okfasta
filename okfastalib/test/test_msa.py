import math
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
    res = column_stats('aaaabbbb--')
    assert res.nvals == 8
    assert res.prop_gaps == 0.2
    assert res.entropy == approx(0.6931472)
    assert res.consensus_val == "a"
    assert res.consensus_prop == 0.5

def test_column_stats_allgaps():
    res = column_stats('----------')
    assert res.nvals == 0
    assert res.prop_gaps == 1.0
    assert res.entropy == 0.0
    assert res.consensus_val == "-"
    assert res.consensus_prop == 1.0

def test_column_stats_nogaps():
    res = column_stats('ccggttaacc')
    assert res.nvals == 10
    assert res.prop_gaps == 0.0
    assert res.entropy == approx(1.332179)
    assert res.consensus_val == "c"
    assert res.consensus_prop == 0.4
