import argparse 

from tblastn_wrapper import tblastn_wrapper

def parse_args():
    parser = argparse.ArgumentParser(
        description="Find ultraconserved elements across multiple genomes"
    )

    parser.add_argument(
        "sequences", 
        metavar="SEQ", 
        type=str, 
        help="the sequence that you want to search"
    )

    parser.add_argument(
        "-out", 
        metavar="OUT", 
        type=str, 
        help="the output location"
    )

    return parser.parse_known_args()

def main():
    args, unknown = parse_args()
    tblastn_wrapper(args, unknown)

if __name__ == "__main__":
    main()
