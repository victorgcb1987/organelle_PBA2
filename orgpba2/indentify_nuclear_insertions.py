from argparse import ArgumentParser
from pathlib import Path

from src.config import OUTPUT_FOLDERS as out_dir
from src.minimap2 import run_minimap2_for_insertions
from src.nuclear_insertions import (get_reads_alignments_info, calculate_reads_query_coverage, 
                                    select_reads_by_coverage, filter_reads_by_name, write_aligned_reads_into_fasta,
                                    get_insertions_positions, exclude_reads_with_less_coverage)


def parse_arguments():
    desc = "Detect organelle insertions from long reads. It needs a organelle assembly and nuclear assembly"
    parser = ArgumentParser(description=desc)

    help_nuclear = "(Required) Reference genome assembly"
    parser.add_argument("--nuclear_assembly", 
                        "-n", type=str,
                        help=help_nuclear,
                        required=True)
    
    help_organelle = "(Required) Organelle genome assembly"
    parser.add_argument("--organelle_assembly", 
                        "-a", type=str,
                        help=help_organelle,
                        required=True)
    
    help_reads = "(Required) Sequences file"
    parser.add_argument("--sequences", 
                        "-s", type=str,
                        help=help_reads,
                        required=True)
    
    help_output = "(Required) Output dir"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    
    help_sequence_tech = "(Required) Sequencing technology used, supported modes: pacbio, pacbio-hifi, nanopore. By default Pacbio"
    parser.add_argument("--technology", "-St",
                        required=True, type=str, 
                        help=help_sequence_tech)
    
    help_threads = "(Optional) Number of threads. 1 thread by default"
    parser.add_argument("--threads", "-t",
                        type=int, help=help_threads,
                        default=1)
    help_coverage_cutoff = "(Optional) coverage cutoff of a read aligned against an organelle to identify potential insertion. By default 0.85)"
    parser.add_argument("--coverage_cutoff", "-c",
                        type=float, help=help_coverage_cutoff,
                        default=0.85)
    help_organelle_comparison = "(Optional) exclude reads with better alignments to this genome. "
    parser.add_argument("--exclude", 
                        "-e", type=str,
                        help=help_organelle_comparison,
                        default="")
    return parser
    

def get_options():
    options = parse_arguments().parse_args()
    nuclear_assembly_fpath = Path(options.nuclear_assembly)
    organelle_assembly_fpath = Path(options.organelle_assembly)
    sequences_fpath = Path(options.sequences)
    output_dir = Path(options.out)
    num_threads = options.threads
    sequence_technology = options.technology
    coverage_cutoff = options.coverage_cutoff
    if options.exclude:
        exclude = Path(options.exclude)
    else:
        exclude = False

    
    return {'nuclear_assembly': nuclear_assembly_fpath,
            'sequences': sequences_fpath, 
            'out_dir': output_dir,
            'number_threads': num_threads,
            'sequence_technology': sequence_technology,
            'organelle_assembly': organelle_assembly_fpath,
            'coverage_cutoff': coverage_cutoff,
            'exclude_assembly': exclude}


def main():
    arguments = get_options()
    if not arguments["out_dir"].exists():
        arguments["out_dir"].mkdir(parents=True)

    organelle_alignments = arguments["out_dir"] / out_dir["minimap2"] /  "mappings_against_organelle.paf"
    arguments["alignment_fpath"] = organelle_alignments
    if not organelle_alignments.exists():
        organelle_alignments = run_minimap2_for_insertions(arguments, assembly="organelle")["output_file"]
    
    if arguments["exclude_assembly"]:
        exclude_alignments =  arguments["out_dir"] / out_dir["minimap2"] /  "mappings_against_organelle_to_exclude.paf"
        arguments["alignment_fpath"] = exclude_alignments
        if not exclude_alignments.exists():
            exclude_alignments = run_minimap2_for_insertions(arguments, assembly="exclude")["output_file"]

    with open(organelle_alignments) as alignments_fhand:
        organelle_alignments_info = get_reads_alignments_info(alignments_fhand)
        if arguments["exclude_assembly"]:
            with open(exclude_alignments) as alignments_fhand:
                exclude_alignments_info = get_reads_alignments_info(alignments_fhand)
                organelle_alignments_info = exclude_reads_with_less_coverage(organelle_alignments_info, exclude_alignments_info)

        organelle_coverages = calculate_reads_query_coverage(organelle_alignments_info)
        reads_under_coverage = select_reads_by_coverage(organelle_coverages, coverage_cutoff=arguments["coverage_cutoff"], mode="under")
        filtered_reads = filter_reads_by_name(organelle_alignments_info, reads_under_coverage)
        seqs_check = arguments["out_dir"] / out_dir["seqs_to_check"]
        if not seqs_check.exists():
            seqs_check.mkdir(parents=True)
        seqs_check = seqs_check / "reads_mapped_against_organelle.fasta"
        if not seqs_check.exists():
            write_aligned_reads_into_fasta(arguments,filtered_reads,seqs_check)
        arguments["sequences"] = seqs_check
    nuclear_alignments = arguments["out_dir"] / out_dir["minimap2"] /  "mappings_against_nuclear.paf"
    if not nuclear_alignments.exists():
        arguments["alignment_fpath"] = nuclear_alignments
        nuclear_alignments = run_minimap2_for_insertions(arguments, assembly="nuclear")["output_file"]
    with open(nuclear_alignments) as alignments_fhand:
        nuclear_alignments_info = get_reads_alignments_info(alignments_fhand)
        reads_over_coverage = select_reads_by_coverage(organelle_coverages, coverage_cutoff=0.9, mode="over")
        for read in reads_over_coverage:
            print(read)
     
    insertions_reads = get_insertions_positions(nuclear_alignments_info,organelle_alignments_info)

        

    #nuclear_alignments = run_minimap2_for_insertions(arguments, organelle=False)


if __name__ == "__main__":
    main()