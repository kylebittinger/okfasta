import argparse
import signal
import sys

from .fasta import parse_fasta, write_fasta
from .util import parse_seq_ids, filter_seq_ids, get_seq_lengths


def filterids_subcommand(args):
    seq_ids = set(parse_seq_ids(args.idsfile))
    seqs = parse_fasta(args.input)
    filtered_seqs = filter_seq_ids(seqs, seq_ids)
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
    subparsers = main_parser.add_subparsers(help='subcommands')

    filterids_parser = subparsers.add_parser(
        "filterids", parents=[common_parser],
        help='Filter by sequence ID')
    filterids_parser.add_argument(
        "idsfile", type=argparse.FileType('r'),
        help="File containing sequence IDs, one per line",
    )
    filterids_parser.set_defaults(func=filterids_subcommand)

    length_parser = subparsers.add_parser(
        "seqlength", parents=[common_parser],
        help='Return sequence lengths in TSV format')
    length_parser.set_defaults(func=length_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None:
        args.input = sys.stdin
    if args.output is None:
        args.output = sys.stdout
    args.func(args)

