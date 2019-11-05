from okfastalib.parse import *

def test_parse_regions_extra_fields():
    f = ["abc 5 7 whatever else"]
    assert list(parse_regions(f)) == [("abc", 5, 7)]

def test_parse_regions():
    f = ["seq1 2 5", "seq3 1 3", "seq1 4 7"]
    regions = [("seq1", 2, 5), ("seq3", 1, 3), ("seq1", 4, 7)]
    assert list(parse_regions(f)) == regions

def test_parse_seq_ids():
    f = ["Id1\n", "\tId2|345 678  \n"]
    assert list(parse_seq_ids(f)) == ["Id1", "Id2|345"]
