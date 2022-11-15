from io import StringIO

def parse_fasta(f):
    f = iter(f)
    line = next(f)
    line = line.strip()
    assert line.startswith(">")
    desc = line[1:]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            seq = StringIO()
        else:
            # Whitespace is not meaningful
            line = line.replace(" ", "")
            seq.write(line)
    yield desc, seq.getvalue()


def write_fasta(f, seqs):
    for desc, seq in seqs:
        f.write(">{0}\n{1}\n".format(desc, seq))
