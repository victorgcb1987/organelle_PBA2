from subprocess import run

from src.dependencies import get_executables
from src.config import OUTPUT_FOLDERS as out_dir
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FILENAMES as out_fname
from src.seqs import get_seq_length
from src.utils import file_exists



def run_subsmapling(options, overwrite=False):
    #Runs subsamplig using filtlong
    #By default, it subsamples by quality (reads sharing kmers with reference genome)
    mapped_sequences = options["seqs_input"]
    orgpba2_run_dir = options["out_dir"]
    desired_coverage = options["desired_coverage"]
    reference_length = get_seq_length(options["reference_input"])
    #Generate number of nucleotides needed for reaching coverage needed
    num_nucleotides = reference_length * desired_coverage
    out_fpath = orgpba2_run_dir / out_dir["subsampled_seqs"] / out_fname["subsampled_seqs"].format(desired_coverage)
    if file_exists(out_fpath) and not overwrite:
        msg = "output {} already exists, skipping"
        msg += " subsampling mapping step"
        results = {"output_file": out_fpath,
                   "return_code": 0,
                   "log_messages": msg.format(str(out_fpath))}
    else:
        filtlong_executable = get_executables(exec_reqs["filtlong"])
        cmd = [filtlong_executable]
        cmd.append("-t {}".format(num_nucleotides))
        cmd.append(str(mapped_sequences))

        if not options["filtlong"]:
            #It adds maximum weight to quality value in order to filter reads for subsampling
            filtlong_options = "-a {} --mean_q_weight 10".format(str(options["reference_input"]))
            cmd.append(filtlong_options)

        cmd.append("| gzip > {}".format(str(out_fpath)))

        filtlong_run = run(" ".join(cmd), shell=True, 
                        capture_output=True)
        results = {"output_file": out_fpath,
                   "return_code": filtlong_run.returncode,
                   "log_messages": filtlong_run.stderr.decode()}
    return results


    
    
