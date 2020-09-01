import re
import sys
import multiprocessing
import tempfile
import argparse
import os
import subprocess

def tblastn_wrapper(args: argparse.Namespace, unknown):
    inp = args.sequences # needs to be modified

    dir_name = tempfile.mkdtemp()

    extension = inp[inp.find('.'):]  # e.g. .fa or .fq
    with open(inp, "r") as f: # reading line by line, failsafe for very large programs
        for line in f:
            while (re.search('>', line) != None):
                name = line
                name = name.rstrip()
                organism_name = name[name.find('>')+2:] #+ extension 
                organism_name = organism_name.replace(" ", "_")
                to_write = tempfile.NamedTemporaryFile(prefix=organism_name, suffix=extension, dir=dir_name, mode="w", delete=False)
                to_write.write(line)
                line = next(f)
                while (re.search('>', line) == None and line != None):
                    to_write.write(line)
                    try:
                        line = next(f)
                    except:
                        break
                to_write.close()

    processes = []

    for file in os.listdir(dir_name):
        file = dir_name + "/" + file 
        cmd = " ".join(unknown)
        command_to_run = ["tblastn", file, cmd]
        p = multiprocessing.Process(target = run_command, args=command_to_run)
        process.append(p)

    for process in processes:
        process.join()

    for file in os.listdir(dir_name):
        os.remove(file)

def run_command(command):
    output = subprocess.run(command, cwd=dir_name, capture_output=True, check=True)

    if (args.out == None):
        print(output)
    else: 
        with open(args.out, "w") as f:
            f.write(output)