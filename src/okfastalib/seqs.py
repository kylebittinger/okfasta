import collections
import itertools
import re


def replace_chars(seqs, replacements):
    for desc, seq in seqs:
        for x, y in replacements:
            seq = seq.replace(x, y)
        yield desc, seq

def replace_seq_ids(seqs, new_seqids):
    for desc, seq in seqs:
        toks = re.split("(\\s)", desc, maxsplit=1)
        old_seq_id = toks[0]
        new_seq_id = new_seqids.get(old_seq_id)
        if new_seq_id is None:
            yield (desc, seq)
        else:
            rest = toks[1:]
            new_desc = new_seq_id + ''.join(rest)
            yield (new_desc, seq)

def get_seq_id(desc):
    return desc.split()[0]

def get_kmers(seqs, k=8):
    for desc, seq in seqs:
        n_kmers = len(seq) - k + 1
        if n_kmers < 1:
            continue
        seq_id = get_seq_id(desc)
        for i in range(n_kmers):
            yield (seq_id, i + 1, seq[i:(i + k)])

def extract_regions(regions, seqs):
    regions_table = collections.defaultdict(list)
    for seq_id, start_pos, end_pos in regions:
        regions_table[seq_id].append((start_pos, end_pos))

    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        for start_pos, end_pos in regions_table[seq_id]:
            extract_id = "{0}__{1}_{2}".format(seq_id, start_pos, end_pos)
            start_idx = start_pos - 1
            end_idx = end_pos
            extract_seq = seq[start_idx:end_idx]
            yield extract_id, extract_seq

def filter_seq_ids(seqs, seq_ids, remove=False):
    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        if remove:
            if seq_id not in seq_ids:
                yield desc, seq
        else:
            if seq_id in seq_ids:
                yield desc, seq

def get_seq_lengths(seqs):
    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        yield seq_id, len(seq)

def search_desc(seqs, regex_str):
    for seq_id, seq in seqs:
        if re.search(regex_str, seq_id):
            yield seq_id, seq

def search_seqs(seqs, query, search_revcomp=False):
    queryset = deambiguate(query)
    if search_revcomp:
        queryset = queryset + [reverse_complement(q) for q in queryset]
    for seq_id, seq in seqs:
        for qseq in queryset:
            if qseq in seq:
                yield seq_id, seq
                continue

def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)

def deambiguate(seq):
    nt_choices = [AMBIGUOUS_BASES[x] for x in seq]
    return ["".join(c) for c in itertools.product(*nt_choices)]

AMBIGUOUS_BASES = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
    }

COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "S": "S",
    "W": "W",
    "H": "D",
    "B": "V",
    "V": "B",
    "D": "H",
    "N": "N",
}
