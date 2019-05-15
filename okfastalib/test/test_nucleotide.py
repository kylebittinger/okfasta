from okfastalib.nucleotide import *

def test_reverse_complement():
    assert reverse_complement("ACGGTAT") == "ATACCGT"

def test_reverse_complement_n():
    assert reverse_complement("AGNN") == "NNCT"

def test_reverse_complement_ambiguous():
    assert reverse_complement("RNSG") == "CSNY"
