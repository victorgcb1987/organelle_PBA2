import gzip
import shutil

from subprocess import run

from src.config import (MAGIC_NUMS_COMPRESSED, 
                        OUTPUT_FOLDERS, 
                        OUTPUT_FILENAMES)

def sequence_kind(path):
    #Checks if sequence file is in fasta or fastq format
    if file_is_compressed(path):
       cat_command = ["zcat"]
    else:
        cat_command = ["cat"]
    cat_command += [str(path), "|", "head", "-n 4"]
    run_cat_command = run("\t".join(cat_command), 
                          shell=True, capture_output=True)
    
    #fastq files headers starts with "@"
    #fasta file headers starts with (">")
    if run_cat_command.stdout.startswith(b'@'):
        return "fastq"
    elif run_cat_command.stdout.startswith(b'>'):
        return "fasta"
    else:
        msg = "{} doesn't seems to be a valid"
        msg += " fasta/fastq file"
        raise RuntimeError(msg.format(path.name))

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
        if not filepath.is_dir() and filepath.stat().st_size > 0:
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

def is_step_done(program, output_dir):
    step_fdir = output_dir / OUTPUT_FOLDERS[program]
    step_fname = step_fdir / OUTPUT_FILENAMES[program]
    if folder_exists(step_fdir) and file_exists(step_fname):
        return True
    else:
        return False

def _compress_file(uncompressed_fpath, compressed_fpath):
    with open(uncompressed_fpath, 'rb') as f_in:
        with gzip.open(compressed_fpath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            f_out.flush()

def compress_file(uncompressed_fpath, output_dir, overwrite=False):
    compressed_file_dir = output_dir / OUTPUT_FOLDERS["mapped_seqs"]
    if not folder_exists(compressed_file_dir):
        compressed_file_dir.mkdir(parents=True, exist_ok=True)
    if sequence_kind(uncompressed_fpath) == "fasta": 
        compressed_file_path = compressed_file_dir / "02_mapped_seqs.fasta.gz"
    if sequence_kind(uncompressed_fpath) == "fastq":
        compressed_file_path = compressed_file_dir / "02_mapped_seqs.fq.gz"
    if file_exists(compressed_file_path) and not overwrite:
        log_messages = "File already exists"
    else:
        _compress_file(uncompressed_fpath, compressed_file_path)
        log_messages = "File compressed"
    return {"output_files": compressed_file_path,
            "log_messages": log_messages}

def check_results(program, results):
    if results["return_code"] == 0:
        print("{} successfully run".format(program))
    else:
        print("{} failed".format(program))
        raise RuntimeError(results["log_messages"])


def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))









    


