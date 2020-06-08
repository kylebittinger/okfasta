import os.path
import tempfile

from okfastalib.command import (
    okfasta_main, msa_ok_main,
)
from okfastalib.fasta import parse_fasta

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

small_fp = os.path.join(data_dir, "small.fasta")
small_aligned_fp = os.path.join(data_dir, "small_aligned.fasta")

rumino_fp = os.path.join(data_dir, "ruminococcus.fasta")
rumino_aligned_fp = os.path.join(data_dir, "ruminococcus_aligned.fasta")


def run_msa_ok(argv, input_fp):
    output_file = tempfile.NamedTemporaryFile(mode='w+t')
    argv.extend(["--input", input_fp, "--output", output_file.name])
    msa_ok_main(argv)
    output_file.seek(0)
    res = output_file.readlines()
    return res

def test_colstats():
    output = run_msa_ok(["colstats"], rumino_aligned_fp)

    assert output[0].startswith("column_position\t")
    assert output[1] == "1\t1\t0.90\t0.0000\tT\t1.00\n"

def test_selectcol_subcommand():
    column_file = tempfile.NamedTemporaryFile(mode="w+t")
    column_file.write("1\n6\n8\n")
    column_file.seek(0)

    output = run_msa_ok(["selectcol", column_file.name], small_aligned_fp)

    assert list(parse_fasta(output)) == [('a|b 42', '-CA'), ('c|2.1 d', 'GCG')]
