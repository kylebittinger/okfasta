import argparse
import random
import signal
import sys

from .seqs import (
    filter_seq_ids, get_seq_lengths, search_seqs, extract_regions,
    search_desc, get_kmers, replace_seq_ids, replace_chars, reverse_complement,
)
from .msa import MSA
from .io import (
    parse_fasta, write_fasta, parse_seq_ids, parse_regions, parse_column_idxs,
    parse_new_ids,
    )

def normalize_subcommand(args):
    seqs = parse_fasta(args.input)
    write_fasta(args.output, seqs)

def replacechars_subcommand(args):
    if args.replace is None:
        replacements = []
    else:
        replacements = [(x, y) for x, y in args.replace]
    if args.remove is not None:
        for x in args.remove:
            replacements.append((x, ''))
    seqs = parse_fasta(args.input)
    replaced_seqs = replace_chars(seqs, replacements)
    write_fasta(args.output, replaced_seqs)

def replaceids_subcommand(args):
    seqs = parse_fasta(args.input)
    new_ids = dict(parse_new_ids(args.newidsfile))
    relabeled_seqs = replace_seq_ids(seqs, new_ids)
    write_fasta(args.output, relabeled_seqs)

def randomseqs_subcommand(args):
    seqs = list(parse_fasta(args.input))
    if args.n > len(seqs):
        args.n = len(seqs)
    selected_seqs = random.sample(seqs, args.n)
    write_fasta(args.output, selected_seqs)

def kmers_subcommand(args):
    seqs = parse_fasta(args.input)
    for seq_id, pos, kmer in get_kmers(seqs, args.k):
        args.output.write("{0}\t{1}\t{2}\n".format(seq_id, pos, kmer))

def extract_subcommand(args):
    seq_regions = parse_regions(args.regionfile)
    seqs = parse_fasta(args.input)
    extracted_seqs = extract_regions(seq_regions, seqs)
    write_fasta(args.output, extracted_seqs)

def revcomp_subcommand(args):
    seqs = parse_fasta(args.input)
    rseqs = ((desc, reverse_complement(seq)) for desc, seq in seqs)
    write_fasta(args.output, rseqs)

def selectcol_subcommand(args):
    seqs = parse_fasta(args.input)
    column_idxs = parse_column_idxs(args.columnfile)
    msa = MSA.from_seqs(seqs)
    msa.filter_by_index(column_idxs, remove=args.remove_columns)
    write_fasta(args.output, msa.seqs)

def colstats_subcommand(args):
    seqs = parse_fasta(args.input)
    msa = MSA.from_seqs(seqs)
    args.output.write("\t".join(msa.column_stats_header))
    args.output.write("\n")
    for stats_result in msa.column_stats():
        stats_values = stats_result.values()
        args.output.write(msa.column_stats_fmt.format(*stats_values))
        args.output.write("\n")

def filterids_subcommand(args):
    seq_ids = set(parse_seq_ids(args.idsfile))
    seqs = parse_fasta(args.input)
    filtered_seqs = filter_seq_ids(seqs, seq_ids, remove=args.remove_ids)
    write_fasta(args.output, filtered_seqs)

def searchdesc_subcommand(args):
    seqs = parse_fasta(args.input)
    filtered_seqs = search_desc(seqs, args.regex)
    write_fasta(args.output, filtered_seqs)

def searchseq_subcommand(args):
    seqs = parse_fasta(args.input)
    filtered_seqs = search_seqs(
        seqs, args.query, search_revcomp=args.search_revcomp)
    write_fasta(args.output, filtered_seqs)

def length_subcommand(args):
    seqs = parse_fasta(args.input)
    for seq_id, seq_len in get_seq_lengths(seqs):
        args.output.write("{0}\t{1}\n".format(seq_id, seq_len))

fasta_io_parser = argparse.ArgumentParser(add_help=False)
fasta_io_parser.add_argument(
    "--input", type=argparse.FileType('r'), default=sys.stdin,
    help="Input FASTA file (default: stdin)",
)
fasta_io_parser.add_argument(
    "--output", type=argparse.FileType('w'), default=sys.stdout,
    help="Output file (default: stdout)",
)

def okfasta_main(argv=None):
    # Ignore SIG_PIPE and don't throw exceptions on it
    # newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='Subcommands')

    extract_parser = subparsers.add_parser(
        "extract", parents=[fasta_io_parser],
        help='Extract sequence regions')
    extract_parser.add_argument(
        "regionfile", type=argparse.FileType('r'),
        help=(
            "File containing sequence ID, start position, and stop "
            "position for each region to extract."))
    extract_parser.set_defaults(func=extract_subcommand)

    filterids_parser = subparsers.add_parser(
        "filterids", parents=[fasta_io_parser],
        help='Filter by sequence ID')
    filterids_parser.add_argument(
        "idsfile", type=argparse.FileType('r'),
        help="File containing sequence IDs, one per line")
    filterids_parser.add_argument(
        "--remove-ids", action="store_true",
        help="Remove, rather than keep, IDs in list")
    filterids_parser.set_defaults(func=filterids_subcommand)

    kmers_parser = subparsers.add_parser(
        "kmers", parents=[fasta_io_parser],
        help='Write k-mers in TSV format')
    kmers_parser.add_argument(
        "--k", type=int, default=8,
        help="K-mer size (default: %(default)s)")
    kmers_parser.set_defaults(func=kmers_subcommand)

    length_parser = subparsers.add_parser(
        "length", parents=[fasta_io_parser],
        help='Write sequence lengths in TSV format')
    length_parser.set_defaults(func=length_subcommand)

    normalize_parser = subparsers.add_parser(
        "normalize", parents=[fasta_io_parser],
        help='Rewrite FASTA file in standard format')
    normalize_parser.set_defaults(func=normalize_subcommand)

    randomseqs_parser = subparsers.add_parser(
        "randomseqs", parents=[fasta_io_parser],
        help='Select random sequences')
    randomseqs_parser.add_argument(
        "--n", type=int, default=100,
        help="Number of sequences (default: %(default)s)")
    randomseqs_parser.set_defaults(func=randomseqs_subcommand)

    replacechars_subparser = subparsers.add_parser(
        "replacechars", parents=[fasta_io_parser],
        help="Replace characters in the sequences")
    replacechars_subparser.add_argument(
        "--replace", type=str, nargs=2, action="append",
        help="Characters to replace")
    replacechars_subparser.add_argument(
        "--remove", type=str, action="append",
        help="Characters to remove")
    replacechars_subparser.set_defaults(func=replacechars_subcommand)

    replaceids_subparser = subparsers.add_parser(
        "replaceids", parents=[fasta_io_parser],
        help="Replace sequence IDs")
    replaceids_subparser.add_argument(
        "newidsfile", type=argparse.FileType('r'),
        help=(
            "File containing existing sequence ID and replacement "
            "sequence ID, one pair per line, separated by whitespace. "
            "Existing sequence IDs not in the file are left as they are."))
    replaceids_subparser.set_defaults(func=replaceids_subcommand)

    revcomp_parser = subparsers.add_parser(
        "revcomp", parents=[fasta_io_parser],
        help='Reverse complement sequences')
    revcomp_parser.set_defaults(func=revcomp_subcommand)

    searchdesc_parser = subparsers.add_parser(
        "searchdesc", parents=[fasta_io_parser],
        help='Find sequences where description matches pattern')
    searchdesc_parser.add_argument(
        "regex",
        help="Regular expression for searching descriptions")
    searchdesc_parser.set_defaults(func=searchdesc_subcommand)

    searchseq_parser = subparsers.add_parser(
        "searchseq", parents=[fasta_io_parser],
        help='Find sequences that match a query sequence exactly')
    searchseq_parser.add_argument(
        "query",
        help="Query sequence")
    searchseq_parser.add_argument(
        "--search-revcomp", action="store_true",
        help="Search for the query and its reverse complement",
    )
    searchseq_parser.set_defaults(func=searchseq_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None: # pragma: no cover
        args.input = sys.stdin
    if args.output is None: # pragma: no cover
        args.output = sys.stdout
    args.func(args)


def msa_ok_main(argv=None):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='Subcommands')

    selectcol_parser = subparsers.add_parser(
        "selectcol", parents=[fasta_io_parser],
        help='Select columns by position')
    selectcol_parser.add_argument(
        "columnfile", type=argparse.FileType('r'),
        help="File containing column numbers, one per line",
    )
    selectcol_parser.add_argument(
        "--remove-columns", action="store_true",
        help="Remove, rather than keep, columns in list",
    )
    selectcol_parser.set_defaults(func=selectcol_subcommand)

    colstats_parser = subparsers.add_parser(
        "colstats", parents=[fasta_io_parser],
        help='Compute summary statistics')
    colstats_parser.set_defaults(func=colstats_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None: # pragma: no cover
        args.input = sys.stdin
    if args.output is None: # pragma: no cover
        args.output = sys.stdout
    args.func(args)

