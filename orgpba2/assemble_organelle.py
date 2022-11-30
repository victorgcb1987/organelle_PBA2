
import argparse
import sys
import json

from pathlib import Path
from uuid import uuid1

from Bio import SeqIO

from src.utils import file_exists, folder_exists

from src.canu import run_canu
from src.config  import EXECUTABLES_REQUIREMENTS, OUTPUT_FILENAMES
from src.config import OUTPUT_FOLDERS as out_dir
from src.heteroplasmy import (find_blocks_breakpoints, create_haplotype_breakpoints_sequences, 
                              write_haplotypes_breakpoints, calculate_heteroplasmy_ratios, 
                              write_heteroplasmy_results)
from src.subsampling import run_subsmapling
from src.dependencies import (check_executable_in_user_envs, 
                              check_is_executable_is_in_PATH)
from src.minimap2 import run_minimap2, run_minimap2_for_heteroplasmy
from src.polishing import run_racon
from src.seqs import (check_if_assembly_is_complete, get_seqs_id_from_paf_file,
                      write_seqs_from_seqs_id,
                      get_seq_length, find_circularity,
                      find_origin, remove_circularity_redundancy,
                      reconstruct_assembly_from_origin,
                      run_blast, get_reads_total_nucleotides, count_seqs,
                      write_circular_region)
from src.utils import (parse_executable_options, 
                       compress_file, check_results, chunkstring)

#Generating program options
def parse_arguments():
    desc = "Assemble organelles from long reads and a reference genome"
    parser = argparse.ArgumentParser(description=desc)
    help_reference = "(Required) Reference genome."
    parser.add_argument("--reference_sequence", 
                        "-r", type=str,
                        help=help_reference,
                        required=True)
    help_sequences = "(Required) Long reads to assemble."
    parser.add_argument("--sequences", "-s",
                        type=str, help=help_sequences,
                        required=True)
    help_output = "(Required) Output dir"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    help_threads = "(Optional) Number of threads. 1 thread by default"
    parser.add_argument("--threads", "-t",
                        type=int, help=help_threads,
                        default=1)
    help_minimap2 = "(Optional) Minimap2 additional arguments. By default excludes secondary alignments"
    parser.add_argument("--minimap2", "-m",
                        type=str, help=help_minimap2,
                        default="secondary=no")
    help_force_mapping = "(Optional) Force minimap2 mapping if it's already done. By default skips this step"
    parser.add_argument("--force_mapping", "-Fm", action='store_true',
                        help=help_force_mapping)
    help_force_extract_mapped_reads = "(Optional) Force extract mapped reads if it's already done. By default skips this step"
    parser.add_argument("--force_extract_mapped_reads", "-Fe", action='store_true',
                        help=help_force_extract_mapped_reads)
    help_sequence_tech = "(Optional) Sequencing technology used, supported modes: pacbio, pacbio-hifi, nanopore. By default Pacbio"
    parser.add_argument("--technology", "-St", default="pacbio", 
                        help=help_sequence_tech)
    help_force_assembly = "(Optional) Force assembly if it is already done. By default skips this step"
    parser.add_argument("--force_assembly", "-fa", action='store_true',
                        help=help_force_assembly)
    help_canu_options = "(Optional) User canu options. None by default"
    parser.add_argument("--canu_options", 
                        type=str, default="",
                        help=help_canu_options)
    help_subsample = "(Optional) Subsample number of mapped reads to desired coverage using filtlong with quality priority set to 10. By default is 0 (subsampling skipped)"
    parser.add_argument("--subsample_coverage", "-sc",
                        type=int,
                        default=0,
                        help=help_subsample)
    help_subsampling = "(Optional) Force subsampling if it is already done. By default skips this step"
    parser.add_argument("--force_subsampling", "-fs", action='store_true',
                        help=help_subsampling)
    help_filtlong = "(Optional) Filtlong user arguments for subsampling. By default quality is given maximum priority"
    parser.add_argument("--filtlong_options",
                        type=str,
                        default="",
                        help=help_filtlong)

    help_polishing = "(Optional) Number of polishing iterations after assembly and before redundancy is removed. 0 by default"
    parser.add_argument("--polishing_iterations",
                        type=int, default=0, help=help_polishing)
    help_racon = "(Optional) Racon polishing additional options. By default is none"
    parser.add_argument("--racon", type=str, default="", help=help_racon)

    help_heteroplasmy = "(Optional) Calculate heteroplasmy if assembly is complete. By default is False"
    parser.add_argument("--heteroplasmy", action='store_true', help=help_heteroplasmy)

    help_metadata = "(optional) get metadata from json file in order to fill fields for results file, if not specified all metadata fields will be set to unidentified"
    parser.add_argument("--metadata", type=Path, default="", help=help_metadata)




    return parser

#Parse and return values given to options when running this program
def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    reference_fpath = Path(options.reference_sequence)
    sequences_fpath = Path(options.sequences)
    output_dir = Path(options.out)
    num_threads = options.threads
    minimap2_options = parse_executable_options(options.minimap2)
    force_mapping = options.force_mapping
    force_extract_mapped_reads = options.force_extract_mapped_reads
    sequence_technology = options.technology
    force_assembly = options.force_assembly
    desired_coverage = options.subsample_coverage
    filtlong_options = options.filtlong_options
    force_subsampling = options.force_subsampling
    num_polishing_iterations = options.polishing_iterations
    calculate_heteroplasmy = options.heteroplasmy
    if options.canu_options:
        canu_options = parse_executable_options(options.canu_options)
    else:
        canu_options = ""
    if options.racon:
        racon_options = parse_executable_options(options.racon)
    else:
        racon_options = ""
    
    
    if options.metadata.is_file():
        metadata_fhand = open(options.metadata)
        metadata = metadata_fhand.read()
        metadata = json.loads(metadata)
        species = metadata.get("species", "NotFound")
        accession = metadata.get("accession", "NotFound")
        sra = metadata.get("SRA", "NotFound")
    else:
        species = "unidentified"
        accession = "unidentified"
        sra = "unidentified"

    return {'reference_input': reference_fpath,
            'sequences': sequences_fpath, 
            'out_dir': output_dir,
            'number_threads': num_threads,
            'minimap2_options': minimap2_options,
            'force_mapping': force_mapping,
            'force_extract_mapped_reads': force_extract_mapped_reads,
            'sequence_technology': sequence_technology,
            'force_flye': force_assembly,
            'desired_coverage': desired_coverage,
            'filtlong': filtlong_options,
            'force_subsampling': force_subsampling,
            'polishing_iterations': num_polishing_iterations,
            'canu_additional_options': canu_options,
            'racon_additional_options': racon_options,
            'calculate_heteroplasmy': calculate_heteroplasmy,
            'species': species,
            'accession': accession,
            'sra': sra}

def main():
    options = get_options()
    output_dir = options["out_dir"]
    reference_fpath = options["reference_input"]

    stats = []
    
    

    #First check if executable requirements are met
    msg = "Check if program dependencies are met:"
    for program, user_env in EXECUTABLES_REQUIREMENTS.items():
        if not check_executable_in_user_envs(user_env, program=program):
            if not check_is_executable_is_in_PATH(user_env["executable"]):
                msg = "{} not found".format(program)
                raise RuntimeError(msg)

    #Get length of the referency assembly used for coverage purpouses
    options["genome_size"] = get_seq_length(reference_fpath)

    #Create output dir
    if not folder_exists(output_dir):
        output_dir.mkdir(parents=True, exist_ok=True)

    #Create log file
    log_number = uuid1()
    log_fhand = open(output_dir / "orgpba2.{}.log".format(log_number), "w")

    command = " ".join(sys.argv)
    msg = "command used: {}\n".format(command)
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    stats.append(msg)

    reads_dataset_number_of_reads = count_seqs(options["sequences"])
    msg = "Total number of reads in dataset: {}\n".format(reads_dataset_number_of_reads)
    stats.append(msg)
    print(msg)
    log_fhand.write(msg)
    reads_dataset_number_of_nucleotides = get_reads_total_nucleotides(options["sequences"])
    msg = "Total number of nucleotides in dataset: {} Gbs\n".format(reads_dataset_number_of_nucleotides)
    stats.append(msg)
    print(msg)
    log_fhand.write(msg)

    #Mapping sequences to reference assembly
    msg = "Step 1: mapping reads against reference sequence with minimap2\n"
    log_fhand.write(msg)
    log_fhand.flush()
    minimap2_results = run_minimap2(options, overwrite=options["force_mapping"])
    check_results(msg, minimap2_results)
    log_fhand.write(minimap2_results["log_messages"])
    log_fhand.flush()
    minimap2_output = minimap2_results["output_file"]

    #Extracting mapped reads
    compressed_output = output_dir / out_dir["mapped_seqs"] / OUTPUT_FILENAMES["mapped_seqs"]
    if not file_exists(compressed_output):
        msg = "Step 2: extracting mapped reads with seqtk\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()
        mapped_seqs_id = get_seqs_id_from_paf_file(minimap2_output)
        #Reads are stored in a temporary file, then compressed
        temp_file = output_dir / "temp"
        mapped_reads_results = write_seqs_from_seqs_id(mapped_seqs_id,  options["sequences"], temp_file, overwrite=options["force_extract_mapped_reads"])
        log_fhand.write(mapped_reads_results["log_messages"])
        log_fhand.flush()
        check_results(msg, mapped_reads_results)
        msg = "Step 2: compressing mapped reads\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()
        compress_file_results = compress_file(temp_file, output_dir, overwrite=options["force_extract_mapped_reads"])
        temp_file.unlink()
        #Downstream steps will use all the mapped reads unless subsampling was set
        options["seqs_input"] = compress_file_results["output_files"]
    else:
        msg = "Step 2: extracting mapped reads with seqtk already done, skipping\n"
        options["seqs_input"] = compressed_output
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()
    
    total_number_of_mapped_reads = count_seqs(options["seqs_input"])
    msg = "Total number of mapped reads: {}\n".format(total_number_of_mapped_reads)
    print(msg)
    stats.append(msg)
    log_fhand.write(msg)
    log_fhand.flush()

    total_number_of_mapped_nucleotides = get_reads_total_nucleotides(options["seqs_input"])
    msg = "Total number of mapped nucleotides: {} Gbs\n".format(total_number_of_mapped_nucleotides)
    print(msg)
    stats.append(msg)
    log_fhand.write(msg)
    log_fhand.flush()

    
    #A subsample will be created from mapped reads if a coverage value was set
    if options["desired_coverage"] > 0:
            msg = "Step 2b: Subsampling to desired coverage of {}\n".format(str(options["desired_coverage"]))
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()
            subsampling_results = run_subsmapling(options, overwrite=options["force_subsampling"])
            check_results(msg, subsampling_results)
            #Subsample will replace all mapped reads for donwstream steps
            
            options["seqs_input"] = subsampling_results["output_file"]
            total_number_of_mapped_reads = count_seqs(options["seqs_input"])
            msg = "Total number of mapped reads, subsampled to coverage {}: {}\n".format(options["desired_coverage"], 
                                                                                         total_number_of_mapped_reads)
            print(msg)
            stats.append(msg)
            log_fhand.write(msg)
            log_fhand.flush()

            total_number_of_mapped_nucleotides = get_reads_total_nucleotides(options["seqs_input"])
            msg = "Total number of mapped nucleotides, subsampled to coverage {}: {} Gbs\n".format(options["desired_coverage"],
                                                                                                   total_number_of_mapped_nucleotides)
            print(msg)
            stats.append(msg)
            log_fhand.write(msg)
            log_fhand.flush()

    # Assembly perfomance with canu
    msg = "Step 3: asssemble contigs with canu\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    canu_results = run_canu(options)
    log_fhand.write(canu_results["log_messages"])
    log_fhand.flush()
    check_results(msg, canu_results)
    options["assembly_fpath"] = canu_results["output_files"]
    options["assembly_fpath"] =  options["assembly_fpath"] / "03_assembled.contigs.fasta"

    # Polishing with racon
    if options["polishing_iterations"] > 0:
        msg = "Step 4: Running_polishing_step, {} iterations".format(options["polishing_iterations"])
        log_fhand.write(msg)
        log_fhand.flush()
        racon_results = run_racon(options)
        log_fhand.write(racon_results["log_messages"])
        log_fhand.flush()
        check_results(msg, racon_results)
        options["assembly_fpath"] = racon_results["output_file"]

    total_number_of_contigs = count_seqs(options["assembly_fpath"])
    msg = "Total number of assembled contigs: {}\n".format(total_number_of_contigs)
    print(msg)
    stats.append(msg)
    log_fhand.write(msg)
    log_fhand.flush()



    #Get the largest contig from assembled contigs. If it's equal or larger than reference assembly, assembly will be considered complete
    #If assembly is complete, circularity will be removed and heteroplasmy calculated
    #If not, script stops here
    assembly_is_incomplete, largest_contig = check_if_assembly_is_complete(options["assembly_fpath"], options["genome_size"], margin=10)
    if assembly_is_incomplete:
        msg = "WARNING: Assembly fragmented! Stopping now\n"
        log_fhand.write(msg)
        stats.append("Assembly fragmented!")
        log_fhand.flush()
    else:
        # Find reference sequence starting point in our assembly 
        no_redundancy_fdir = output_dir / out_dir["no_redundant"]
        no_redundancy_fpath = no_redundancy_fdir / OUTPUT_FILENAMES["no_redundant"]
        redundant_seq_fpath = no_redundancy_fdir / "04_redundant_sequence.fasta"
        if not folder_exists(no_redundancy_fdir):
            no_redundancy_fdir.mkdir(parents=True, exist_ok=True)
        if not file_exists(no_redundancy_fpath):
            with open(redundant_seq_fpath, "w") as redundant_fhand: 
                SeqIO.write(largest_contig, redundant_fhand, "fasta")
        blast_arguments = {"kind": "blastn",
                           "task": "megablast",
                           "outmft": "-outfmt '6 std qlen slen'"}
        origin_blastn = run_blast(redundant_seq_fpath, reference_fpath,
                              options, blast_arguments=blast_arguments)
        origin = find_origin(origin_blastn["output_file"], margin=10)
        circularity = find_circularity(redundant_seq_fpath, options, overlap_length=100)
        if circularity:
            circularity_length = abs(circularity["overlap_start"][0]-circularity["overlap_start"][1])
            redundant_seq_length = get_seq_length(redundant_seq_fpath)
            msg = "Length of assembly with circularity: {}\n".format(redundant_seq_length)
            msg += "Circularity found at {}\n".format(circularity)
            msg += "Length of the circularity: {}\n".format(circularity_length)
            log_fhand.write(msg)
            stats.append(msg)
            log_fhand.flush()
            circular_fpath = output_dir / "05_remove_circularity_redundancy" / "04_circular_sequence.fasta"
            write_circular_region(redundant_seq_fpath, circularity, circular_fpath)
            msg = "Circular sequence write in {}\n".format(circular_fpath)
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()

        # Remove redundance by circularituy
        no_redundant_seq = remove_circularity_redundancy(redundant_seq_fpath, circularity)
        reconstructed_sequence = reconstruct_assembly_from_origin(no_redundant_seq, origin)
        #Write sequence
        no_redundancy_fhand = open(no_redundancy_fpath, "w")
        no_redundancy_fhand.write(">no_redundant_assembly\n")
        sequence_chunks = chunkstring(reconstructed_sequence, 70)
        no_redundancy_fhand.write("\n".join(sequence_chunks))
        no_redundancy_fhand.close()
        reference_seq_length = get_seq_length(reference_fpath)
        assembly_length = get_seq_length(no_redundancy_fpath)
        msg = "Length of the reference genome used for this assembly: {}\nLength of the assembly created: {}\n".format(str(reference_seq_length), str(assembly_length))
        print(msg)
        stats.append(msg)
        log_fhand.write(msg)
        log_fhand.flush()
    results = output_dir / "orgpba2.{}.stats".format(log_number)
    with open(results, "w") as stats_fhand:
        for stat in stats:
            stats_fhand.write(stat)
            stats_fhand.flush()

    # Calculate heteroplasmy
    if options["calculate_heteroplasmy"]:
        breakpoints, colinear = find_blocks_breakpoints(no_redundancy_fpath, options)
        print(breakpoints)
        sequence = open(no_redundancy_fpath)
        sequence.readline()
        sequence = "".join([line.rstrip() for line in sequence])
        plastid_parts = create_haplotype_breakpoints_sequences(sequence, breakpoints, colinear=colinear)
        options["haplotypes_fpath"] = write_haplotypes_breakpoints(options, plastid_parts, overwrite=True)
        options["mapped_reads_against_haplotypes_fpath"] = options["out_dir"] / out_dir["haplotypes"] / "mapped_reads_against_haplotypes.paf"
        options["seqs_input"] = compressed_output
        results = run_minimap2_for_heteroplasmy(options)
        alignment_results = open(results["output_file"])
        heteroplasmy_fhand = open(options["out_dir"] / "heteroplasmy_results.txt")
        heteroplasmy_results = calculate_heteroplasmy_ratios(alignment_results, breakpoints)
        write_heteroplasmy_results(heteroplasmy_results, heteroplasmy_fhand)


if __name__ == "__main__":
    main()