import argparse
from os import cpu_count

from tblastn_wrapper.tblastn_wrapper import tblastn_wrapper


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find ultraconserved elements across multiple genomes"
    )

    parser.add_argument(
        "-query",
        metavar="QUERY",
        type=str,
        required=True,
        help="the sequence that you want to search",
    )

    parser.add_argument(
        "-out", metavar="OUT", type=str, default="output.fa", help="the output location"
    )

    parser.add_argument(
        "-threads",
        metavar="THREADS",
        type=int,
        help="The number of threads to use",
        default=cpu_count(),
    )

    return parser.parse_known_args()


def main():
    args, unknown = parse_args()

    tblastn_wrapper(args, unknown)


if __name__ == "__main__":
    main()
