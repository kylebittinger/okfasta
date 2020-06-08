import os.path
import tempfile

from okfastalib.command import (
    okfasta_main, msa_ok_main,
)

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
rumino_fp = os.path.join(data_dir, "ruminococcus.fasta")
rumino_aligned_fp = os.path.join(data_dir, "ruminococcus_aligned.fasta")

def test_colstats():
    output_file = tempfile.NamedTemporaryFile(mode='w+t')
    argv = [
        "colstats",
        "--input", rumino_aligned_fp,
        "--output", output_file.name]
    msa_ok_main(argv)
    output_file.seek(0)
    print(output_file.read())
    output_file.seek(0)    

    header_line = output_file.readline()
    assert header_line.startswith("column_position\t")

    col1_line = output_file.readline()
    assert col1_line == "1\t3\t"
