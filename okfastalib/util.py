def get_seq_id(desc):
    return desc.split()[0]

def parse_seq_ids(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        seq_id = line.split()[0]
        yield seq_id

def filter_seq_ids(seqs, seq_ids):
    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        if seq_id in seq_ids:
            yield desc, seq

def get_seq_lengths(seqs):
    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        yield seq_id, len(seq)

