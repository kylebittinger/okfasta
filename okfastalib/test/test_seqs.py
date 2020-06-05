from okfastalib.seqs import *

def test_get_seq_id():
    assert get_seq_id("AB.C|13260 DEF ghijk") == "AB.C|13260"

def test_select_columns():
    seqs = [("hg", "ACGTT-"), ("s2", "--GTAC")]
    assert list(select_columns(seqs, [1, 3, 4, 5])) == [
        ("hg", "AGTT"), ("s2", "-GTA")]

def test_tabulate_positions():
    seqs = [("kj dfus", "ACG"), ("fh:9", "-T")]
    rows = tabulate_positions(seqs)
    assert list(rows) == [
        ("kj", 1, "A"), ("kj", 2, "C"), ("kj", 3, "G"),
        ("fh:9", 1, "-"), ("fh:9", 2, "T"),
    ]

def test_extract_regions():
    regions = [("seq1", 2, 5)]
    seqs = [("seq1", "ACGACTA")]
    assert list(extract_regions(regions, seqs)) == [("seq1__2_5", "CGAC")]

def test_search_desc():
    seqs = [("abcde", "GCTTG"), ("PTR1.1 genA", "CTCTCG")]
    assert list(search_desc(seqs, r"\wcd")) == [seqs[0]]
    assert list(search_desc(seqs, r"^[^b]+$")) == [seqs[1]]

def test_search_seqs():
    seqs = [("n", "CGTTAC"), ("m", "CTGGTGTCA")]
    assert list(search_seqs(seqs, "GTTA")) == [seqs[0]]

def test_search_seqs_revcomp():
    seqs = [("n", "CGTTAC"), ("m", "CTGGTGTCA")]
    assert list(search_seqs(seqs, "ACAC", True)) == [seqs[1]]

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
