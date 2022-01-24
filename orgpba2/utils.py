from pathlib import Path

from orgpba2.config import MAGIC_NUMS_COMPRESSED

def parse_args(args):
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

def file_exists(filepath, filename):
    filepath = Path(filepath) / filename
    return filepath.exists()

def file_is_compressed(path):
    max_len = max(len(x) for x in MAGIC_NUMS_COMPRESSED)
    with open(path, 'rb') as fhand:
        file_start = fhand.read(max_len)
        fhand.close()
    for magic_num in MAGIC_NUMS_COMPRESSED:
        if file_start.startswith(magic_num):
            return True
    return False




    


