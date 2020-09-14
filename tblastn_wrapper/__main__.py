import argparse 
import subprocess
import sys
from os import cpu_count
from tblastn_wrapper.wrapper_function import wrapper


def parse_args():

    path_command = subprocess.run('which tblastn', 
                        check=True, 
                        shell=True, 
                        capture_output=True
                        )

    tblastn_path = path_command.stdout.decode("utf-8").partition('\n')[0]
    version_command = tblastn_path + " -version"
    version_res = subprocess.run(version_command, shell=True, capture_output=True)

    parser = argparse.ArgumentParser(
        add_help=False, usage="tblastn_wrapper [commands from tblastn]"
    )

    parser.add_argument(
        "-help",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-h",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-version",
        action="version",
        help=argparse.SUPPRESS,
        version=argparse._(version_res.stdout.decode("utf-8"))
    )

    parser.add_argument(
        "-query", 
        metavar="QUERY", 
        type=str, 
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-out", 
        metavar="OUT", 
        help=argparse.SUPPRESS,
        type=str, 
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
    path_command = subprocess.run('which tblastn', 
                                    shell=True, 
                                    capture_output=True
                                )

    if (path_command.returncode != 0):
        raise Exception("tblastn is not installed or initialised")

    args, unknown = parse_args()

    tblastn_path = path_command.stdout.decode("utf-8").partition('\n')[0]
    help_command = tblastn_path + " -help"
    h_command = tblastn_path + " -h"

    help_res = subprocess.run(help_command, shell=True, capture_output=True)
    h_res = subprocess.run(h_command, shell=True, capture_output=True)
    
    if args.h:
        if (args.help):
            print(help_res.stdout.decode("utf-8"))
        else:
            print(h_res.stdout.decode("utf-8"))
        exit(0)
    elif args.help:
        print(help_res.stdout.decode("utf-8"))
        exit(0)

    if "-db" not in unknown or args.query is None:
        sys.stderr.write("BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified\nPlease refer to the BLAST+ user manual.")
        exit(path_command.returncode)

    wrapper(args, unknown)

if __name__ == "__main__":
    main()