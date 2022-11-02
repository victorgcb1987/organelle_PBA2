from ctypes import alignment
from src.config import OUTPUT_FILENAMES as out_fnames
from src.config import OUTPUT_FOLDERS as out_dirs
from src.seqs import complement, get_seq_length, reverse_complement, run_blast
from src.utils import file_exists, folder_exists


def find_blocks_breakpoints(sequence, options):
    #It finds breakpints between single copy fragments and Inverted Repeats in the reference genome/Assembly
    #It also detects if Inverted repeats are colinear
    blast_arguments = {"kind": "blastn",
                        "task": "megablast",
                        "outmft": "-outfmt '6 std qlen slen'" }
    sequence_length = get_seq_length(sequence)
    blastn_results = run_blast(sequence, sequence, 
                               options, blast_arguments)
    blast_results_output = blastn_results["output_file"]
    with open(blast_results_output) as blast_ouput:
        breakpoint = {}
        for line in blast_ouput:
            line = line.split()
            aln_length = int(line[3])
            reference_start = int(line[6])
            reference_end = int(line[7])
            query_start = int(line[8])
            query_end = int(line[9])
            if aln_length == sequence_length:
                continue
            else:
                if not breakpoint:
                    breakpoint = {"ref_start": reference_start,
                                  "ref_end": reference_end,
                                  "query_start": query_start,
                                  "query_end": query_end,
                                  }
                #Ira and irb have inverted orientation
                elif query_start  == breakpoint["ref_end"] and query_end == breakpoint["ref_start"]:
                    return {"LSC_IRa": reference_start, "IRa_SSC": reference_end,  
                            "SSC_IRb": query_end, "IRb_LSC": query_start}, False
                #Ira and irb are colinear          
                elif query_start == breakpoint["ref_start"] and query_end == breakpoint["ref_end"]:
                    return {"LSC_IRa": reference_start, "IRa_SSC": reference_end,  
                            "SSC_IRb": query_start, "IRb_LSC": query_end}, True

def create_haplotype_breakpoints_sequences(sequence, breakpoints, colinear=False):
    initial_position = 0
    parts = {}
    for value, position in breakpoints.items():
        #converting value to 0 level representation
        position -= 1
        #get the the nearest higher breakpoint to the one we are creating
        higher_breakpoints = []
        for next_value, next_position in breakpoints.items():
            if position < next_position-1:
                higher_breakpoints.append({"name": next_value, 
                                           "position": next_position -1})

        #Get different orientations of parts
        print(higher_breakpoints)
        first_part_name, _ = value.split("_")
        print(value, initial_position, position)
        first_part_regular = sequence[initial_position:position]
        first_part_reversed = sequence[initial_position:position][::-1]
        first_part_complement = complement(sequence[initial_position:position])
        first_part_revcomp = reverse_complement(sequence[initial_position:position])   
        initial_position = position
        print(first_part_name)
        parts[first_part_name] = {"forward": first_part_regular, 
                                  "reverse": first_part_reversed, 
                                  "complement": first_part_complement,
                                  "revcomp": first_part_revcomp}
        if first_part_name == "LSC":
            parts[first_part_name].pop('reverse')
            parts[first_part_name].pop('revcomp')
        if first_part_name == "IRa" and colinear:
            parts[first_part_name]["colinear"] = first_part_regular
    return parts


def write_haplotypes_breakpoints(options, parts, overwrite=False):
    out_dir = options["out_dir"] / out_dirs["haplotypes"]
    out_fpath = out_dir / out_fnames["haplotypes"]
    if not folder_exists(out_dir):
        out_dir.mkdir(parents=True, exist_ok=True)
    if file_exists(out_fpath) and not overwrite:
        return out_fpath
    #Reconstructs all possible haplotypes, with LSC as starting point
    out_fhand = open(out_fpath, "w")
    first_parts = parts["LSC"]
    second_parts = parts["IRa"]
    third_parts = parts["SSC"]
    for orientation1, sequence1 in first_parts.items():
        for orientation2, sequence2, in second_parts.items():
            for orientation3, sequence3 in third_parts.items():
                if orientation2 == "forward":
                    orientation4 = "revcomp"
                    sequence4 = second_parts[orientation4]
                elif orientation2 == "complement":
                    orientation4 = "reverse"
                    sequence4 = second_parts[orientation4]
                elif orientation2 == "reverse":
                    orientation4 = "complement"
                    sequence4 = second_parts[orientation4]
                elif orientation2 == "revcomp":
                    sequence4 = second_parts[orientation4]
                elif orientation4 == "colinear":
                    sequence4 = second_parts["forward"]
                    orientation2 = "forward"
                    orientation4 = "forward"
                sequence = ">{}_{}-{}_{}-{}_{}-{}_{}\n".format("LSC", orientation1, 
                                                               "IR", orientation2,
                                                               "SSC", orientation3,
                                                               "IR", orientation4)
                haplotype = sequence1 +  sequence2 + sequence3 + sequence4
                duplicated_haplotype = sequence + haplotype + haplotype + "\n"

                out_fhand.write(duplicated_haplotype)
                out_fhand.flush()
    out_fhand.close()
    return out_fpath


def calculate_heteroplasmy_ratios(alignment_results_fhand, breakpoints, min_distance_covered=1000):
    reads_found_in_haplotypes = {}
    repetitive_regions_positions = [(breakpoints["LSC_IRa"], breakpoints["IRa_SSC"]),  
                                    (breakpoints["SSC_IRb"], breakpoints["IRb_LSC"]),
                                    (breakpoints["LSC_IRa"] + breakpoints["IRb_LSC"], breakpoints["IRa_SSC"] + breakpoints["IRb_LSC"]),  
                                    (breakpoints["SSC_IRb"] + breakpoints["IRb_LSC"], breakpoints["IRb_LSC"] + breakpoints["IRb_LSC"])]
    for read_aligned in alignment_results_fhand:
        read_aligned = read_aligned.split()
        read_id = read_aligned[0]
        haplotype = read_aligned[5]
        alignment_start = int(read_aligned[7])
        alignment_end = int(read_aligned[8])
        for repetitive_region_position in repetitive_regions_positions:
            sequence_start = repetitive_region_position[0] - min_distance_covered
            sequence_end = repetitive_region_position[1] + min_distance_covered
            if alignment_start <= sequence_start and alignment_end >= sequence_end:
                if haplotype not in reads_found_in_haplotypes:
                    reads_found_in_haplotypes[haplotype] = [read_id]
                else:
                    reads_found_in_haplotypes[haplotype].append(read_id)
    return reads_found_in_haplotypes


def write_heteroplasmy_results(heteroplasmy_results, results_fhand):
    for haplotype, reads in heteroplasmy_results.items():
        results_fhand.write("Number of reads supporting haplotype {}: {}\n".format(haplotype, len(reads)))
        results_fhand.write("Reads supporting haplotype {}:\n".format(haplotype))
        results_fhand.write("\n".join(reads)+"\n\n")
        results_fhand.flush()
    results_fhand.close()
 