from subprocess import run

from orgpba2.utils import file_is_compressed

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

def count_seqs(path):
    #count number of seqs in file
    #by counting number of headers
    if file_is_compressed(path):
        grep_command = ["zgrep"]
    else:
        grep_command = ["grep"]

    if sequence_kind(path) == "fastq":
        grep_command += ["\"@\"", str(path), 
                         "|", "wc -l"]
        print(grep_command)
    elif sequence_kind(path) == "fasta":
        grep_command += ["\">\"", str(path),
                         "|", "wc -l"]
    run_grep_command = run("\t".join(grep_command),
                           shell=True, capture_output=True)
    return int(run_grep_command.stdout.decode())

    


