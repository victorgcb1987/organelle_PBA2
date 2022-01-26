import argparse

from pathlib import Path

from orgpba2.config import MAGIC_NUMS_COMPRESSED

def parse_executable_options(args):
    #parse arguments provided by the user at orgpba2 init
    #arguments must have the following syntax : arg1=value1,arg2=value3
    parsed_args = {}
    args = args.split(",")
    for arg in args:
        arg = arg.split("=")
        if len(arg) != 2:
            raise RuntimeError("arguments must be separated by a comma")
        else:
            parsed_args[arg[0]] = arg[1]
    return parsed_args

def file_exists(filepath):
    if filepath.exists():
        if not filepath.is_dir():
            return True
    return False

def folder_exists(filepath):
    if filepath.exists():
        if filepath.is_dir():
            return True
    return False
        

def file_is_compressed(path):
    max_len = max(len(x) for x in MAGIC_NUMS_COMPRESSED)
    with open(path, 'rb') as fhand:
        file_start = fhand.read(max_len)
        fhand.close()
    for magic_num in MAGIC_NUMS_COMPRESSED:
        if file_start.startswith(magic_num):
            return True
    return False

def create_arguments():
    prog= "Organelle_PBA2"
    desc = "An organelle assembler using"
    desc += " using pacbio long reads"
    parser = argparse.ArgumentParser(prog=prog,
                                     description=desc)
    parser.add_argument("-r", type=str)
    parser.add_argument("-i", type=str)
    parser.add_argument("-b", type=str)
    parser.add_argument("-t", type=int)
    return parser


    


