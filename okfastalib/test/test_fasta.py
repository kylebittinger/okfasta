from okfastalib.fasta import *

def test_parse_fasta():
    f = [
        ">ab c", "GG VU", "LG",
        ">bu\tter", "CGTA"
    ]
    seqs = parse_fasta(f, trim_desc=True)
    assert list(seqs) == [("ab", "GGVTLG"), ("bu", "CGTA")]
