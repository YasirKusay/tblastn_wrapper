import re
import multiprocessing
import tempfile
import argparse
import os
import subprocess

def wrapper(args: argparse.Namespace, extra_args):
    query_filename = args.query

    with tempfile.TemporaryDirectory() as dir_name:
        subqueries = split_query(query_filename, dir_name)

        run_queries(subqueries, args.out, dir_name, args.threads, extra_args)

def split_query(query_filename, working_dir):
    subqueries = []

    # Get the file extension (e.g. .fa or .fq)
    _, extension = os.path.splitext(query_filename)

    with open(query_filename, "r") as f:
        for line in f:
            # Search until reaching a FASTA header
            while re.search('>', line) is not None:
                # Extract the query header name
                name = line.rstrip()

                letters = re.search(r"[^ >]", line, re.I)
                start = letters.start()

                to_write = tempfile.NamedTemporaryFile(
                    suffix=extension, 
                    dir=working_dir, 
                    mode="w", 
                    delete=False,
                )

                subqueries.append(to_write.name)

                to_write.write(line)

                # Copy query lines until reaching the next header or blank line
                line = next(f)
                while re.search('>', line) is None and line is not None:
                    to_write.write(line)
                    try:
                        line = next(f)
                    except:
                        break
                to_write.close()

    return subqueries

def run_queries(query_filenames, output_filename, working_dir, threads, extra_args):
    
    query_commands = []

    cmd = " ".join(extra_args)

    # below is in case the user decides to alias the tblastn_wrapper into tblastn
    # in the case of that, we must know the path of the tblastn command
    # done in shell, as locate is very fast

    # command below allows tblastn_wrapper to be
    #  compatible with any version of tblastn
    path_command = subprocess.run('which tblastn', check=True, shell=True, capture_output=True)
    tblastn_path = path_command.stdout.decode("utf-8").partition('\n')[0]

    for filename in query_filenames:
        command = tblastn_path + " -query " + filename + " " + cmd
        query_commands.append(command)

    with multiprocessing.Pool(processes=threads) as pool:
        query_results = pool.map(run_worker, query_commands)

        num_queries = 0

        closing = "" # stores the info, such as matrix used, gap penalties, etc

        if output_filename is not None:
            with open(output_filename, "wb") as f:
                for res in query_results:
                    if (num_queries == 0):
                        num_queries += 1
                        to_print, closing = initial_line_process(res)
                        f.write(to_print)
                    else:
                        f.write(process_lines(res))
                
                f.write(closing)
        else:
            for res in query_results:
                if (num_queries == 0):
                    num_queries += 1
                    to_print, closing = initial_line_process(res)
                    print(to_print.decode("utf-8"))
                else:
                    print(process_lines(res).decode("utf-8"))
            
            print(closing.decode("utf-8"))

    for filename in query_filenames:
        os.remove(filename)

def initial_line_process(lines):
    str_array = lines.decode("utf-8").splitlines()
    str_array.reverse()
    end = 0

    for line in str_array:
        if re.search('^[ ]*Database:', line) is not None:
            break
        end += 1

    str_array.reverse()
    to_return = str_array[:len(str_array)-end-1]
    closing = str_array[len(str_array)-end-1:]

    return "\n".join(to_return).encode('utf-8'), "\n".join(closing).encode('utf-8')


# getting rid of the introductory/conclusion lines that are part of the 
# output of tblastn

def process_lines(lines):
    str_array = lines.decode("utf-8").splitlines()
    start = 0
    for line in str_array:
        if re.search('^Query=', line) is not None:
            break
        start += 1

    # reversing the list is more efficient than iterating through the list
    end = 0 
    str_array.reverse()
    for line in str_array:
        if re.search('^[ ]*Database:', line) is not None:
            break
        end += 1
    
    str_array.reverse()
    to_return = str_array[start:len(str_array)-end-1]
    return "\n".join(to_return).encode('utf-8')


def run_worker(command):
    res = subprocess.run(command, check=True, shell=True, capture_output=True)

    return res.stdout