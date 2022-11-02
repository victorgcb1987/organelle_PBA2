import argparse
from pathlib import Path

from src.config import EXECUTABLES_REQUIREMENTS
from src.dependencies import (check_executable_in_user_envs, 
                              check_is_executable_is_in_PATH)
from src.heteroplasmy import (find_blocks_breakpoints, create_haplotype_breakpoints_sequences, 
                              write_haplotypes_breakpoints, calculate_heteroplasmy_ratios, 
                              write_heteroplasmy_results)
from src.minimap2 import run_minimap2_for_heteroplasmy
from src.utils import folder_exists

def parse_arguments():
    desc = "Calculate heteroplasmy ratios found in a long-read dataset for a given plastid genome"
    parser = argparse.ArgumentParser(description=desc)
    help_reference = "Mandatory plastid assembly in fasta format"
    parser.add_argument("--assembly" , 
                        "-a", type=str,
                        help=help_reference,
                        required=True)
    help_sequences = "Mandatory long reads to assemble"
    parser.add_argument("--sequences", "-s",
                        type=str, help=help_sequences,
                        required=True)
    help_output = "Output dir"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    help_threads = "Number of threads"
    
    parser.add_argument("--threads", "-t",
                        type=int, help=help_threads,
                        default=1)

    help_sequence_tech = "Sequencing technology used, supported modes: pacbio, pacbio-hifi, nanopore"
    parser.add_argument("--technology", "-St", default="pacbio", type=str, 
                        help=help_sequence_tech)
    return parser





def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    reference_fpath = Path(options.assembly)
    sequences_fpath = Path(options.sequences)
    output_dir = Path(options.out)
    num_threads = options.threads
    sequence_technology = options.technology
    

    return {'reference_input': reference_fpath,
            'seqs_input': sequences_fpath, 
            'out_dir': output_dir,
            'number_threads': num_threads,
            'sequence_technology': sequence_technology,}



def main():
    options = get_options()
    output_dir = options["out_dir"]
    if not folder_exists(output_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
    reference_fpath = options["reference_input"]
    msg = "Check if program dependencies are met:"
    for program, user_env in EXECUTABLES_REQUIREMENTS.items():
        if program in ["minimap2", "blastn", "makeblastdb"]:
            if not check_executable_in_user_envs(user_env, program=program):
                if not check_is_executable_is_in_PATH(user_env["executable"]):
                    msg = "{} not found".format(program)
                    raise RuntimeError(msg)
    breakpoints, colinear = find_blocks_breakpoints(reference_fpath, options)
    print(breakpoints)
    sequence = open(reference_fpath)
    sequence.readline()
    sequence = "".join([line.rstrip() for line in sequence])
    plastid_parts = create_haplotype_breakpoints_sequences(sequence, breakpoints, colinear=colinear)
    options["haplotypes_fpath"] = write_haplotypes_breakpoints(options, plastid_parts, overwrite=True)
    results = run_minimap2_for_heteroplasmy(options)
    alignment_results = open(results["output_file"])
    heteroplasmy_results = calculate_heteroplasmy_ratios(alignment_results, breakpoints)
    results_fhand = open(options["out_dir"] / "heteroplasmy_results.txt", "w")
    write_heteroplasmy_results(heteroplasmy_results, results_fhand)
    





if __name__ == "__main__":
    main()