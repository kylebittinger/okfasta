import collections
import re

from .nucleotide import reverse_complement, deambiguate

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
