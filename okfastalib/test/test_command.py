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
