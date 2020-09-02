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

        run_queries(subqueries, args.out, dir_name, extra_args)


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


def run_queries(query_filenames, output_filename, working_dir, extra_args):
    processes = []

    for file in query_filenames:
        cmd = " ".join(extra_args)
        command_to_run = "tblastn -query " + file + " " + cmd
        p = multiprocessing.Process(
            target=run_command, args=(command_to_run, output_filename, working_dir)
        )
        p.start()
        processes.append(p)

    for process in processes:
        process.join()


def run_command(command, output_filename, working_dir):
    if output_filename is None:
        subprocess.run(command, cwd=working_dir, check=True, shell=True)
    else:
        with open(output_filename, "a") as f:
            subprocess.run(command, stdout=f, cwd=working_dir, check=True, shell=True)
