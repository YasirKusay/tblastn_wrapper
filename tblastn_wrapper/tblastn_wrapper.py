import re
import multiprocessing
import tempfile
import argparse
import os
import subprocess


def tblastn_wrapper(args: argparse.Namespace, extra_args):
    query_filename = args.query

    with tempfile.TemporaryDirectory() as dir_name:
        subqueries = split_query(query_filename, dir_name)

        run_queries(subqueries, args.out, dir_name, args.threads, extra_args)


def split_query(query_filename, working_dir):
    subqueries = []

    # Get the file extension (e.g. .fa or .fq)
    _, extension = os.path.splitext(query_filename)

    with open(query_filename, "r") as f:
        # reading line by line, failsafe for very large queries
        for line in f:
            # Search until reaching a FASTA header
            while re.search(">", line) is not None:
                # Extract the query header name
                name = line.rstrip()

                letters = re.search(r"[^ >]", line, re.I)
                start = letters.start()

                organism_name = name[name.find(">") + start :]  # + extension
                organism_name = organism_name.replace(" ", "_")

                # Split the query into a new temporary file
                to_write = tempfile.NamedTemporaryFile(
                    prefix=organism_name,
                    suffix=extension,
                    dir=working_dir,
                    mode="w",
                    delete=False,
                )

                subqueries.append(to_write.name)

                to_write.write(line)

                # Copy query lines until reaching the next header or blank line
                line = next(f)
                while re.search(">", line) is None and line is not None:
                    to_write.write(line)
                    try:
                        line = next(f)
                    except:
                        break
                to_write.close()

    return subqueries


def run_queries(query_filenames, output_filename, working_dir, threads, extra_args):
    query_commands = [
        ["tblastn", "-query", filename, extra_args] for filename in query_filenames
    ]

    with multiprocessing.Pool(processes=threads) as pool:
        query_results = pool.map(run_worker, query_commands)

        with open(output_filename, "wb") as f:
            for res in query_results:
                f.write(res)


def run_worker(command):
    res = subprocess.run(command, check=True, shell=True, capture_output=True)

    return res.stdout
