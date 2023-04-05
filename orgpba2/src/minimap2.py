from os import environ as env
from subprocess import run

from src.config import OUTPUT_FILENAMES as out_fname
from src.config import OUTPUT_FOLDERS as out_dir
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.dependencies import get_executables
from src.utils import file_exists, folder_exists

def _set_output_filepath(run_minimap2_dir, 
                        is_sam_output=False):
    if is_sam_output:
        out_filepath = run_minimap2_dir / "01_mapping_minimap2.sam"
    else:
        out_filepath = run_minimap2_dir / "01_mapping_minimap2.paf"
    
    return out_filepath
    

def run_minimap2(options, overwrite=False, haplotypes=False):
    #Init required variables
    reference_input = options["reference_input"]
    sequences = options["sequences"]
    additional_arguments = options["minimap2_options"]
    orgpba2_run_dir = options["out_dir"]
    run_minimap2_dir = orgpba2_run_dir / out_dir["minimap2"]
    
    #Check whether minimap2 executable is in user env
    #or path and start minimap's run command adding it
    minimap2_executable = get_executables(exec_reqs["minimap2"])
    cmd = [minimap2_executable]
    if options["sequence_technology"] == "pacbio":
        cmd.append("-x map-pb")
    if options["sequence_technology"] == "nanopore":
        cmd.append("-x map-ont")
    if options["sequence_technology"] == "pacbio-hifi":
        cmd.append("-x map-hifi")

    #minimap2 can output in PAF/SAM format
    #if -a is set it will output unfiltered SAM
    #else it will output aligned reads only in PAF
    is_sam_output = False
    for argument, value in additional_arguments.items():
        #overriding user minimap2 output option
        #for piping between programs purpouses
        if argument == "o":
            continue
        elif "ax" in argument:
            is_sam_output = True
            cmd.append("-{} {}".format(argument, value))
        elif argument == "secondary":
            cmd.append("--secondary {}".format(value))
        
    #If num_threads is set globally when running org_pba2
    #and -t is not set in minimap2 arguments we will set
    #number of threads to global value
    if "number_threads" in options:
        if "-t" not in additional_arguments:
            cmd += ["-t {}".format(str(options["number_threads"]))]
    
    #adding reference sequence and pacbio reads
    cmd += ["{} {}".format(reference_input, sequences)]

    #setting output file path
    out_filepath = _set_output_filepath(run_minimap2_dir, 
                                        is_sam_output)
    if not folder_exists(run_minimap2_dir):
        run_minimap2_dir.mkdir(parents=True, exist_ok=True)
    if file_exists(out_filepath) and not overwrite:
        msg = "output {} already exists, skipping"
        msg += " minimap2 mapping step"
        results = {"output_file": out_filepath,
                   "return_code": 0,
                   "log_messages": msg.format(str(out_filepath))}
        return results
    if is_sam_output:
        cmd += ["|samtools view -F 2052"]
    cmd += ["> {}".format(str(out_filepath))]

    #running minimap2
    minimap2_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)

    #return needed info for checking if program has run 
    #properly, log messages and files generated
    results = {"output_file": out_filepath,
               "return_code": minimap2_run.returncode,
               "log_messages": minimap2_run.stderr.decode()}
    return results

def run_minimap2_for_polishing(options):
    #Minimap run for polishing
    minimap2_executable = get_executables(exec_reqs["minimap2"])
    cmd = [minimap2_executable]
    if options["sequence_technology"] == "pacbio":
        cmd.append("-x map-pb")
    if options["sequence_technology"] == "nanopore":
        cmd.append("-x map-ont")
    if options["sequence_technology"] == "pacbio-hifi":
        cmd.append("-x map-hifi")
    cmd.append(str(options["assembly_to_polish_fpath"]))
    cmd.append(str(options["seqs_input"]))
    cmd.append("-t {}".format(str(options["number_threads"])))
    cmd.append("> {}".format(str(options["mapped_reads_against_assembly_to_polish_fpath"])))
    minimap2_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)
    results = {"output_file": options["mapped_reads_against_assembly_to_polish_fpath"],
               "return_code": minimap2_run.returncode,
               "log_messages": minimap2_run.stderr.decode()}
    return results


def run_minimap2_for_heteroplasmy(options):
    #Minimap run for heteroplasmy ratio calculations
    minimap2_executable = get_executables(exec_reqs["minimap2"])
    cmd = [minimap2_executable]
    output_fpath = options["out_dir"] / out_dir["haplotypes"] / out_fname["mapping_for_haplotypes"]
    if file_exists(output_fpath):
        results = results = {"output_file": output_fpath,
                             "return_code": 0,
                             "log_messages": "File already exists"}
        return results
    if options["sequence_technology"] == "pacbio":
        cmd.append("-x map-pb")
    if options["sequence_technology"] == "nanopore":
        cmd.append("-x map-ont")
    if options["sequence_technology"] == "pacbio-hifi":
        cmd.append("-x map-hifi")
    cmd.append("--secondary=no -L -c")
    cmd.append(str(options["haplotypes_fpath"]))
    cmd.append(str(options["seqs_input"]))
    cmd.append("-t {}".format(str(options["number_threads"])))
    cmd.append("> {}".format(str(output_fpath)))
    minimap2_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)
    results = {"output_file": output_fpath,
               "return_code": minimap2_run.returncode,
               "log_messages": minimap2_run.stderr.decode()}
    return results


def run_minimap2_for_insertions(options, assembly=""):
    #Minimap run for heteroplasmy ratio calculations
    minimap2_executable = get_executables(exec_reqs["minimap2"])
    cmd = [minimap2_executable]
    alns_fdir = output_fpath = options["out_dir"] / out_dir["minimap2"]
    if not alns_fdir.exists():
        alns_fdir.mkdir(parents=True)

    output_fpath = options["alignment_fpath"]
    if file_exists(output_fpath):
        results = results = {"output_file": output_fpath,
                             "return_code": 0,
                             "log_messages": "File already exists"}
        print(results)
        return results
    if options["sequence_technology"] == "pacbio":
        cmd.append("-cx map-pb")
    if options["sequence_technology"] == "nanopore":
        cmd.append("-cx map-ont")
    if options["sequence_technology"] == "pacbio-hifi":
        cmd.append("-cx map-hifi")
    cmd.append("--secondary=no --cs")
    if assembly == "organelle":
        cmd.append(str(options['organelle_assembly']))
    elif assembly == "nuclear":
        cmd.append(str(options['nuclear_assembly']))
    elif assembly == "exclude":
        cmd.append(str(options['exclude_assembly']))
    cmd.append(str(options["sequences"]))
    cmd.append("-t {}".format(str(options["number_threads"])))
    cmd.append("> {}".format(str(output_fpath)))
    minimap2_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)
    results = {"output_file": output_fpath,
               "return_code": minimap2_run.returncode,
               "log_messages": minimap2_run.stderr.decode()}
    print(cmd)
    print(results)
    return results

