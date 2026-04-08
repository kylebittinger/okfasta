from okfastalib.io import *

def test_parse_fasta():
    f = [
        ">ab c", "GG VU", "LG",
        ">bu\tter", "CGTA"
    ]
    seqs = parse_fasta(f)
    assert list(seqs) == [("ab c", "GGVULG"), ("bu\tter", "CGTA")]

def test_parse_column_idxs():
    f = ["3", "   ", " # comment line", "6 with extra stuff"]
    assert list(parse_column_idxs(f)) == [3, 6]

def test_parse_regions_extra_fields():
    f = ["abc 5 7 whatever else"]
    assert list(parse_regions(f)) == [("abc", 5, 7)]

def test_parse_regions():
    f = ["seq1 2 5", "seq3 1 3", " ", "# comment", "seq1 4 7"]
    regions = [("seq1", 2, 5), ("seq3", 1, 3), ("seq1", 4, 7)]
    assert list(parse_regions(f)) == regions

def test_parse_seq_ids():
    f = ["Id1\n", "\tId2|345 678  \n", "   ", "   # a comment", "  id3"]
    assert list(parse_seq_ids(f)) == ["Id1", "Id2|345", "id3"]

def test_parse_new_ids():
    f = ["seq1 jfk", "kj\tjdf\tss", " # comment!", "\t\t\t", "seq3     k"]
    assert list(parse_new_ids(f)) == [
        ("seq1", "jfk"), ("kj", "jdf"), ("seq3", "k")]
