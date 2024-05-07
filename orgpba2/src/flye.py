import shutil

from subprocess import run

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FOLDERS as out_dir
from src.dependencies import get_executables
from src.utils import folder_exists


def run_flye(options, overwrite=False):
    orgpba2_run_dir = options["out_dir"]
    run_flye_dir = orgpba2_run_dir / out_dir["flye"]
    #check if output exists
    if folder_exists(run_flye_dir):
        if overwrite:
            shutil.rmtree(run_flye_dir)
        else:
            msg = "output dir {} already exists, skipping"
            msg += " flye assembling step"
            results = {"output_files": run_flye_dir / "assembly.fasta",
                       "return_code": 0,
                       "log_messages": msg.format(str(run_flye_dir))}
            return results
    #Init required variables
    #additional_arguments = options["flye_options"]
    #Check whether flye executable is in user env
    #or path and start minimap's run command adding it
    flye_executable = get_executables(exec_reqs["flye"])
    cmd = [flye_executable]
    #adding user defined arguments
    if options["sequence_technology"] == "pacbio":
        cmd.append("--pacbio-raw")
    elif options["sequence_technology"] == "pacbio-hifi":
        cmd.append("--pacbio-hifi")
    elif options["sequence_technology"] == "ont":
        cmd.append("--nano-raw")
    else:
        raise RuntimeError("Unsuported sequencing technlogy: {}".format(options["sequence_technology"]))
    cmd.append(str(options["seqs_input"]))
    cmd.append("--genome-size {}".format(options["genome_size"]))
    # if additional_arguments:
    #     for argument, value in additional_arguments.items():
    #         if argument == "min-overlap":
    #             cmd.append("--min-overlap {}".format(value))
    #         if argument == "read_error":
    #             cmd.append("--read_error {}".format(value))
    #         if argument == "scaffold":
    #             cmd.append("--scaffold")
    #         if argument == "i" or argument == "iterations":
    #             cmd.append("--iterations {}".format(value))
    #         if argument == "asm-coverage":
    #             cmd.append("--asm-coverage {}".format(value))
    cmd.append("-t {} --scaffold --iterations 3".format(options["number_threads"]))
    cmd.append("--out-dir {}".format(str(run_flye_dir)))
    #running flye
    print(" ".join(cmd))
    print(cmd)
    flye_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)
    #return needed info for checking if program has run 
    #properly, log messages and files generated
    results = {"output_files": run_flye_dir / "assembly.fasta",
               "return_code": flye_run.returncode,
               "log_messages": flye_run.stderr.decode()}
    return results