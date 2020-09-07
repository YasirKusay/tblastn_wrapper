import argparse 
from os import cpu_count
from tblastn_wrapper.wrapper_function import wrapper

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run multiple tblastn commands at the same time"
    )

    parser.add_argument(
        "-query", 
        metavar="QUERY", 
        type=str, 
        required=True,
        help="the sequence that you want to search"
    )

    parser.add_argument(
        "-out", 
        metavar="OUT", 
        type=str, 
        help="the output location"
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
    wrapper(args, unknown)


if __name__ == "__main__":
    main()
