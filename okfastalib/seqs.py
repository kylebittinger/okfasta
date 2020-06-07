import collections
import re

from .nucleotide import reverse_complement, deambiguate

def get_seq_id(desc):
    return desc.split()[0]

def tabulate_positions(seqs):
    for desc, seq in seqs:
        seq_id = get_seq_id(desc)
        for n, val in enumerate(seq):
            yield (seq_id, n + 1, val)

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
