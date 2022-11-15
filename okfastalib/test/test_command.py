import os.path
import tempfile

from okfastalib.command import (
    okfasta_main, msa_ok_main,
)
from okfastalib.fasta import parse_fasta

small_fasta = """\
>a|b 42
GCAGACGATAC
>c|2.1 d
GCAGCCGGT
"""

tall_fasta = """\
>a
CTT
>b
GAT
>c
TCA
>d
ATG
>e
CAG
>f
ACC
>g
TGG
>h
GCT
>i
TAA
>j
AGC
"""

small_aligned_fasta = """\
>a|b 42
-CAGACGAT-
>c|2.1 d
GC-GCCGG--
"""

def parse_fasta_list(f):
    return list(parse_fasta(f))

def tempfile_containing(contents):
    f = tempfile.NamedTemporaryFile(mode="w+t")
    f.write(contents)
    f.seek(0)
    return f

def run_msa_ok(argv, input_data):
    input_file = tempfile_containing(input_data)
    output_file = tempfile_containing("")
    argv.extend(["--input", input_file.name, "--output", output_file.name])
    msa_ok_main(argv)
    output_file.seek(0)
    res = output_file.readlines()
    return res

def test_colstats():
    output = run_msa_ok(["colstats"], small_aligned_fasta)
    assert output[0].startswith("column_position\t")
    assert output[1] == "1\t1\t0.50\t0.0000\tG\t1.00\n"

def test_selectcol_subcommand():
    column_file = tempfile_containing("1\n6\n8\n")
    output = run_msa_ok(["selectcol", column_file.name], small_aligned_fasta)
    assert list(parse_fasta(output)) == [('a|b 42', '-CA'), ('c|2.1 d', 'GCG')]

def run_okfasta(argv, input_data):
    input_file = tempfile_containing(input_data)
    output_file = tempfile_containing("")
    argv.extend(["--input", input_file.name, "--output", output_file.name])
    okfasta_main(argv)
    output_file.seek(0)
    res = output_file.readlines()
    return res

def test_replacechars_subcommand():
    output = run_okfasta([
        "replacechars",
        "--replace", "A", "B",
        "--replace", "C", "D",
        "--remove", "G",
        ], small_fasta)
    output_seqs = parse_fasta_list(output)
    assert output_seqs == [("a|b 42", "DBBDBTBD"), ("c|2.1 d", "DBDDT")]

def test_replaceids_subcommand():
    newids_file = tempfile_containing("c|2.1\tc-2\naaa ggg\n")
    output = run_okfasta(["replaceids", newids_file.name], small_fasta)
    output_seq_ids = [seq_id for seq_id, seq in parse_fasta_list(output)]
    assert output_seq_ids == ["a|b 42", "c-2 d"]

def test_randomseqs_subcommand():
    output = run_okfasta(["randomseqs", "--n", "3"], tall_fasta)
    input_seqs = parse_fasta_list(tall_fasta.splitlines())
    output_seqs = parse_fasta_list(output)
    assert len(output_seqs) == 3
    for output_seq in output_seqs:
        assert output_seq in input_seqs

def test_randomseqs_subcommand_largen():
    output = run_okfasta(["randomseqs", "--n", "20"], tall_fasta)
    input_seqs = parse_fasta_list(tall_fasta.splitlines())
    output_seqs = parse_fasta_list(output)
    assert list(sorted(output_seqs)) == input_seqs

def test_filterids_subcommand():
    ids_file = tempfile_containing("c|2.1")
    output = run_okfasta(["filterids", ids_file.name], small_fasta)
    assert parse_fasta_list(output) == [("c|2.1 d", "GCAGCCGGT")]

def test_extract_subcommand():
    region_file = tempfile_containing(
        "a|b\t3\t7\n"
        "c|2.1\t1\t2\n")
    output = run_okfasta(["extract", region_file.name], small_fasta)
    assert parse_fasta_list(output) == [
        ("a|b__3_7", "AGACG"),
        ("c|2.1__1_2", "GC"),
    ]

def test_searchdesc_subcommand():
    output = run_okfasta(["searchdesc", "^a"], small_fasta)
    assert parse_fasta_list(output) == [("a|b 42", "GCAGACGATAC")]

def test_seqrchseq_subcommand():
    output = run_okfasta(["search", "AGACGAT"], small_fasta)
    assert parse_fasta_list(output) == [("a|b 42", "GCAGACGATAC")]

def test_length_subcommand():
    output = run_okfasta(["length"], small_fasta)
    assert output == ["a|b\t11\n", "c|2.1\t9\n"]

def test_revcomp_subcommand():
    output = run_okfasta(["revcomp"], small_fasta)
    assert parse_fasta_list(output) == [
        ("a|b 42", "GTATCGTCTGC"),
        ("c|2.1 d", "ACCGGCTGC"),
    ]

def test_kmers_subcommand():
    output = run_okfasta(["kmers", "--k", "8"], small_fasta)
    assert output == [
        "a|b\t1\tGCAGACGA\n",
        "a|b\t2\tCAGACGAT\n",
        "a|b\t3\tAGACGATA\n",
        "a|b\t4\tGACGATAC\n",
        "c|2.1\t1\tGCAGCCGG\n",
        "c|2.1\t2\tCAGCCGGT\n",
    ]
