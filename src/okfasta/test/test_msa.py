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

def test_pairwise_mismatches():
    seqs = [
        ("a b", "GTCC"),
        ("t|k", "G-CA"),
        ("ryt", "AAAA"),
    ]
    assert list(pairwise_mismatches(seqs)) == [
        ("a", "t|k", 1),
        ("a", "ryt", 4),
        ("t|k", "ryt", 2),
    ]
    assert list(pairwise_mismatches(seqs, include_gaps=True)) == [
        ("a", "t|k", 2),
        ("a", "ryt", 4),
        ("t|k", "ryt", 3),
    ]
    assert list(pairwise_mismatches(seqs, percent=True)) == [
        ("a", "t|k", 100 * 1 / 3),
        ("a", "ryt", 100 * 4 / 4),
        ("t|k", "ryt", 100 * 2 / 3),
    ]

def test_mismatches():
    assert mismatches("ABCD", "ABEE") == 2
    assert mismatches("ABCD", "AB-E") == 1
    assert mismatches("ABCD", "AB-E", include_gaps=True) == 2
