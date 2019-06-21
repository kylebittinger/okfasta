from okfastalib.util import *

def test_search_seqs():
    seqs = [("n", "CGTTAC"), ("m", "CTGGTGTCA")]
    assert list(search_seqs(seqs, "GTTA")) == [("n", "CGTTAC")]

def test_search_seqs_revcomp():
    seqs = [("n", "CGTTAC"), ("m", "CTGGTGTCA")]
    assert list(search_seqs(seqs, "ACAC", True)) == [("m", "CTGGTGTCA")]

def test_get_seq_id():
    assert get_seq_id("AB.C|13260 DEF ghijk") == "AB.C|13260"

def test_parse_seq_ids():
    f = ["Id1\n", "\tId2|345 678  \n"]
    assert list(parse_seq_ids(f)) == ["Id1", "Id2|345"]

def test_filter_seq_ids():
    seqs = [("a", "TCGCT"), ("b", "GGCT"), ("cd", "TCGACG")]
    ids = ["cd", "b"]
    # Original order of seqs is retained
    # Only first sequence is missing
    assert list(filter_seq_ids(seqs, ids)) == seqs[1:]

def test_remove_seq_ids():
    seqs = [("a", "TCGCT"), ("b", "GGCT"), ("cd", "TCGACG")]
    ids = ["b"]
    # Original order of seqs is retained
    # Second sequence is missing
    assert list(filter_seq_ids(seqs, ids, remove=True)) == [seqs[0], seqs[2]]

def test_get_seq_lengths():
    seqs = [("ab cde", "GCTCGCT"), ("f|g hij", "GCTCGAGTCA")]
    assert list(get_seq_lengths(seqs)) == [("ab", 7), ("f|g", 10)]
