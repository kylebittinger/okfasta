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

def parse_column_idxs(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        toks = line.split(maxsplit=1)
        column_idx = int(toks[0])
        yield column_idx

def parse_regions(f):
    for n, line in enumerate(f):
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        toks = line.split()
        seq_id = toks[0]
        start_pos = int(toks[1])
        end_pos = int(toks[2])
        yield seq_id, start_pos, end_pos

def parse_seq_ids(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        seq_id = line.split()[0]
        yield seq_id

def parse_new_ids(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        toks = line.split()
        old_seq_id = toks[0]
        new_seq_id = toks[1]
        yield old_seq_id, new_seq_id
