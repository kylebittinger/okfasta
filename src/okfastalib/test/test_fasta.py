from okfastalib.fasta import *

def test_parse_fasta():
    f = [
        ">ab c", "GG VU", "LG",
        ">bu\tter", "CGTA"
    ]
    seqs = parse_fasta(f)
    assert list(seqs) == [("ab c", "GGVULG"), ("bu\tter", "CGTA")]
