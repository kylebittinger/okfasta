import csv

def parse_columns(f):
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

def write_tsv(f, header, rows):
    w = csv.writer(f, lineterminator="\n", delimiter="\t")
    w.writerow(header)
    w.writerows(rows)

def parse_seq_ids(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        seq_id = line.split()[0]
        yield seq_id

