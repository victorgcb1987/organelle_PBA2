import shutil

from subprocess import run

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FOLDERS as out_dir
from src.dependencies import get_executables
from src.utils import folder_exists

def run_canu(options, overwrite=False):
    #Init required variables
    additional_arguments = options["canu_additional_options"]
    orgpba2_run_dir = options["out_dir"]
    run_canu_dir = orgpba2_run_dir / out_dir["canu"]
    #Check whether canu executable is in user env
    #or path and start minimap's run command adding it
    canu_executable = get_executables(exec_reqs["canu"])
    cmd = [canu_executable]
    #adding user defined arguments
    add_prefix = True
    if additional_arguments:
        for argument, value in additional_arguments.items():
            if argument == "assemble":
                cmd.append("-assemble")
            if argument.startswith("-p"):
                add_prefix = False
            if argument == "stopOnLowCoverage":
                cmd.append("{}={}".format(argument, value))
            else:
                cmd.append("{} {}".format(argument, value))
    if add_prefix:
        cmd.append("-p 03_assembled")
    #adding output_dir
    cmd.append("-d {}".format(str(run_canu_dir)))
    cmd.append("-{} {}".format(options["sequence_technology"], options["seqs_input"]))
    cmd.append("genomeSize={}".format(options["genome_size"]))

    #If num_threads is set globally when running org_pba2
    #and -t is not set in minimap2 arguments we will set
    #number of threads to global value
    if "number_threads" in options:
        if "maxThreads" not in additional_arguments:
            cmd += ["maxThreads={}".format(str(options["number_threads"]))]
    if folder_exists(run_canu_dir):
        if overwrite:
            shutil.rmtree(run_canu_dir)
        else:
            msg = "output dir {} already exists, skipping"
            msg += " canu assembling step"
            results = {"output_files": run_canu_dir,
                       "return_code": 0,
                       "log_messages": msg.format(str(run_canu_dir))}
            return results
    

    #running canu
    print(" ".join(cmd))
    canu_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)

    #return needed info for checking if program has run 
    #properly, log messages and files generated
    results = {"output_files": run_canu_dir,
               "return_code": canu_run.returncode,
               "log_messages": canu_run.stderr.decode()}
    return results