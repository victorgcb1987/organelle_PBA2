from subprocess import run


from src.config import OUTPUT_FOLDERS as out_dir
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FILENAMES as out_fname
from src.dependencies import get_executables
from src.minimap2 import run_minimap2_for_polishing
from src.utils import file_exists, folder_exists


def run_racon(options, overwrite=False):
    #It runs racon for polishing step
    #In iteration 1, it will use the original assembly and reads mapped to it
    #In following iterations, it will use the previous polished assembly and reads mapped to it
    racon_executable = get_executables(exec_reqs["racon"])
    assembly_to_polish_fpath = options["assembly_fpath"]
    orgpba2_run_dir = options["out_dir"]
    polishing_dir = orgpba2_run_dir / out_dir["polished_assembly"]
    if not folder_exists(polishing_dir):
        polishing_dir.mkdir(parents=True, exist_ok=True)
    iterations = 0
    while iterations < options["polishing_iterations"]:
        polished_assembly_fpath = polishing_dir / out_fname["polished_assembly"].format(iterations+1)
        polish_mapping_fpath = polishing_dir / out_fname["mapping_for_polishing"].format(iterations+1)
        options["mapped_reads_against_assembly_to_polish_fpath"] = polish_mapping_fpath
        options["assembly_to_polish_fpath"] = assembly_to_polish_fpath
        if not file_exists(polished_assembly_fpath) or overwrite:
            run_minimap2_for_polishing(options)
            sequences_to_map = options["seqs_input"]
            cmd = [racon_executable, str(sequences_to_map), str(polish_mapping_fpath), str(assembly_to_polish_fpath)]
            if options["racon_additional_options"]:
                for additional_option in options["racon_additional_options"]:
                    cmd.append(additional_option)
            cmd.append("> {}".format(polished_assembly_fpath))
            racon_run = run(" ".join(cmd), shell=True, 
                            capture_output=True)
            assembly_to_polish_fpath = polished_assembly_fpath
            return_code = racon_run.returncode
            msg = racon_run.stderr.decode()
        else:
            assembly_to_polish_fpath = polished_assembly_fpath
            return_code = 0
            msg = "polishing step already done"
        iterations += 1
    results = {"output_file": assembly_to_polish_fpath,
               "return_code": return_code,
               "log_messages": msg}
    return results
    



    

