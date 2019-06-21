import argparse
import signal
import sys

from .fasta import parse_fasta, write_fasta
from .util import (
    parse_seq_ids, filter_seq_ids, get_seq_lengths, search_seqs,
)
from .nucleotide import reverse_complement


def revcomp_subcommand(args):
    seqs = parse_fasta(args.input)
    rseqs = ((desc, reverse_complement(seq)) for desc, seq in seqs)
    write_fasta(args.output, rseqs)

def filterids_subcommand(args):
    seq_ids = set(parse_seq_ids(args.idsfile))
    seqs = parse_fasta(args.input)
    filtered_seqs = filter_seq_ids(seqs, seq_ids, remove=args.remove_ids)
    write_fasta(args.output, filtered_seqs)

def search_subcommand(args):
    seqs = parse_fasta(args.input)
    filtered_seqs = search_seqs(
        seqs, args.queryseq, remove=args.search_revcomp)
    write_fasta(args.output, filtered_seqs)

def length_subcommand(args):
    seqs = parse_fasta(args.input)
    for seq_id, seq_len in get_seq_lengths(seqs):
        args.output.write("{0}\t{1}\n".format(seq_id, seq_len))

def main(argv=None):
    # Ignore SIG_PIPE and don't throw exceptions on it
    # newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--input", type=argparse.FileType('r'), default=sys.stdin,
        help="Input FASTA file (default: stdin)",
    )
    common_parser.add_argument(
        "--output", type=argparse.FileType('w'), default=sys.stdout,
        help="Output file (default: stdout)",
    )

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='Subcommands')

    filterids_parser = subparsers.add_parser(
        "filterids", parents=[common_parser],
        help='Filter by sequence ID')
    filterids_parser.add_argument(
        "idsfile", type=argparse.FileType('r'),
        help="File containing sequence IDs, one per line",
    )
    filterids_parser.add_argument(
        "--remove-ids", action="store_true",
        help="Remove, rather than keep, IDs in list",
    )
    filterids_parser.set_defaults(func=filterids_subcommand)

    search_parser = subparsers.add_parser(
        "search", parents=[common_parser],
        help='Find sequences matching query (exact matches)')
    filterids_parser.add_argument(
        "queryseq",
        help="Query sequence",
    )
    filterids_parser.add_argument(
        "--search-revcomp", action="store_true",
        help="Search for the query or its reverse complement",
    )
    filterids_parser.set_defaults(func=search_subcommand)

    length_parser = subparsers.add_parser(
        "seqlength", parents=[common_parser],
        help='Return sequence lengths in TSV format')
    length_parser.set_defaults(func=length_subcommand)

    revcomp_parser = subparsers.add_parser(
        "revcomp", parents=[common_parser],
        help='Reverse complement sequences')
    revcomp_parser.set_defaults(func=revcomp_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None:
        args.input = sys.stdin
    if args.output is None:
        args.output = sys.stdout
    args.func(args)

