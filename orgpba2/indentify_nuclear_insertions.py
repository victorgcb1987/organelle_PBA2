from argparse import ArgumentParser
from pathlib import Path
from statistics import median
from sys import argv

from src.config import OUTPUT_FOLDERS as out_dir
from src.heteroplasmy import find_blocks_breakpoints
from src.minimap2 import run_minimap2_for_insertions
from src.nuclear_insertions import (get_reads_alignments_info, calculate_reads_query_coverage, 
                                    select_reads_by_coverage, filter_reads_by_name, write_seqs_from_seqs_id,
                                    get_insertions_positions, exclude_reads_with_less_coverage, remove_organelle_offset,
                                    group_reads_of_same_insertion)
from src.seqs import (calculate_repeats_regions, concate_reference_genome)


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
    help_organelle_length = "(Required) organelle length"
    parser.add_argument("--length",
                        '-l', type=int,
                        help=help_organelle_length,
                        required=True)
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
    organelle_length = options.length
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
            'exclude_assembly': exclude,
            'organelle_length': organelle_length}


def main():
    arguments = get_options()
    if not arguments["out_dir"].exists():
        arguments["out_dir"].mkdir(parents=True)
    breakpoints, colinear = find_blocks_breakpoints(arguments["organelle_assembly"], arguments)
    repeats = [(breakpoints["LSC_IRa"], breakpoints["IRa_SSC"]), (breakpoints["SSC_IRb"], breakpoints["IRb_LSC"])]
    print(repeats)
    print(colinear)
    arguments["nuclear_assembly"] = concate_reference_genome(arguments["nuclear_assembly"], arguments["out_dir"])
    exit()

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
        organelle_alignments_info = remove_organelle_offset(organelle_alignments_info, arguments['organelle_length'])
        if arguments["exclude_assembly"]:
            with open(exclude_alignments) as alignments_fhand:
                exclude_alignments_info = get_reads_alignments_info(alignments_fhand, exclude_potential_chimeras=False)
                organelle_alignments_info = exclude_reads_with_less_coverage(organelle_alignments_info, exclude_alignments_info)

        organelle_coverages = calculate_reads_query_coverage(organelle_alignments_info)
        reads_under_coverage = select_reads_by_coverage(organelle_coverages, coverage_cutoff=arguments["coverage_cutoff"], mode="under")
        filtered_reads = filter_reads_by_name(organelle_alignments_info, reads_under_coverage)
        seqs_check = arguments["out_dir"] / out_dir["seqs_to_check"]
        if not seqs_check.exists():
            seqs_check.mkdir(parents=True)
        seqs_check = seqs_check / "reads_mapped_against_organelle.fasta.gz"
        if not seqs_check.exists():
            write_seqs_from_seqs_id(filtered_reads,arguments["sequences"],seqs_check)
        arguments["sequences"] = seqs_check
    nuclear_alignments = arguments["out_dir"] / out_dir["minimap2"] /  "mappings_against_nuclear.paf"
    if not nuclear_alignments.exists():
        arguments["alignment_fpath"] = nuclear_alignments
        nuclear_alignments = run_minimap2_for_insertions(arguments, assembly="nuclear")["output_file"]
    with open(nuclear_alignments) as alignments_fhand:
        nuclear_alignments_info = get_reads_alignments_info(alignments_fhand)
        nuclear_coverages = calculate_reads_query_coverage(nuclear_alignments_info)
        reads_over_coverage = select_reads_by_coverage(nuclear_coverages, coverage_cutoff=0.9, mode="over")
        nuclear_alignments_info = {readname: values for readname, values in  nuclear_alignments_info.items() if readname in reads_over_coverage}
    chroms, insertions_reads = get_insertions_positions(nuclear_alignments_info,organelle_alignments_info)
    insertions = group_reads_of_same_insertion(insertions_reads)
    results_fpath = arguments["out_dir"] / "insertions_found.tsv"
    with open(results_fpath, "w") as results_fhand:
        results_fhand.write("#{}".format(" ".join(argv)))
        results_fhand.write("#CHROM_ID\tNUCLEAR_START\tNUCLEAR_END\tORGANELLE_START\tORGANELLE_END\tINSERTION_LENGTH\tNUM_READS\tREADS_ID\n")
        for insertion in insertions:
            chrom = insertion["nuclear"]
            nuclear_start = str(round(median(insertion["insertion_starts"])))
            nuclear_end = str(round(median(insertion["insertion_ends"])))
            organelle_starts = str(round(median(insertion["organelle_starts"])))
            organelle_ends = str(round(median(insertion["organelle_ends"])))
            nupt_length = str(round(median(insertion["organelle_ends"])) - round(median(insertion["organelle_starts"])))
            num_reads = len(round(insertion["readnames"]))
            readnames = ",".join(insertion["readnames"])
            results_fhand.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, nuclear_start, nuclear_end, organelle_starts, organelle_ends, nupt_length, num_reads, readnames))
            results_fhand.flush()



        

    #nuclear_alignments = run_minimap2_for_insertions(arguments, organelle=False)


if __name__ == "__main__":
    main()