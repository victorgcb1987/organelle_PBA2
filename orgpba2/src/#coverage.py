from csv import DictReader

def sort_sequences_by_position(paf_file):
    sorted_alignments = {}
    field_names = ["query_ID", "query_length", "query_start",
                   "query_end", "strand", "subject_ID",
                   "subject_length", "subject_start",
                   "subject_end", "nucleotide_matches", 
                   "aln_length", "mapQ", "1", "2", "3", 
                   "4", "5", "6"]
    with open(paf_file) as paf_fhand:
        for read in DictReader(paf_fhand, fieldnames=field_names, delimiter="\t"):
            for field, value in read.items():
                try:
                    read[field] = int(value)
                except:
                    continue
            if read["nucleotide_matches"] < 1000:
                continue
            if read["query_ID"] in sorted_alignments:
                if read["query_start"] < sorted_alignments[read["query_ID"]]["query_start"]:
                    sorted_alignments[read["query_ID"]] = read
            else:
                 sorted_alignments[read["query_ID"]] = read
    sorted_alignments = [values for values in sorted_alignments.values()]
    sorted_alignments = sorted(sorted_alignments, key=lambda x: x["subject_start"])
    return sorted_alignments

def get_reads_by_coverage(sorted_reads_by_position, reference_length, desired_coverage, fraction_window=0.9):
    position = (0, 100)
    reads_retrieved = []
    coverages = {}
    keep_searching = True
    while keep_searching:
        window_reads = []
        reads = []
        for sorted_read in sorted_reads_by_position:    
            if sorted_read["subject_start"] >= position[0] and sorted_read["subject_start"] <= position[1]:
                reads.append(sorted_read)
        if not reads:
            break
        reads = sorted(reads, key=lambda x: x["subject_end"], reverse=True)
        position = (position[0], int(reads[0]["subject_end"]))
        coverage_to_reach = ((position[1] - position[0])*desired_coverage)
        nucleotides_aligned = 0
        for read in reads:
            aln_start = read["subject_start"]
            if read["subject_end"] > position[1]:
                aln_end = position[1]
            else:
                aln_end = read["subject_end"]
            nucleotides_aligned += abs(aln_end - aln_start)
            window_reads.append(read)
            #sorted_reads_by_position.remove(read)
            if nucleotides_aligned >= coverage_to_reach:
                break
        for read in window_reads:
            sorted_reads_by_position.remove(read)
        if nucleotides_aligned != 0:
            print(position)
            coverage = int(nucleotides_aligned / (position[1] - position[0]))
            total_length_aligned_nucls = 0
            lengths = []
            for win_read in window_reads:
                total_length_aligned_nucls += win_read["aln_length"]
                lengths.append(win_read["aln_length"])
                max_length = 0
                for length in lengths:
                    if length > max_length:
                        max_length = length
            mean_length = int(total_length_aligned_nucls / len(reads))
            coverages[position] = {"coverage": coverage,
                                   "mean_aligned_nucls": mean_length,
                                   "num_reads": len(window_reads),
                                   "max_aligned_length": max_length}
            #offset = int((max_length)*fraction_window)

            position = (int(position[-1]*fraction_window), position[-1]+1000)
        else:
            coverages[position] = {"coverage": 0,
                                   "mean_length": 0,
                                   "num_reads": 0}               
            position = (int(position[-1]*fraction_window), position[-1]+(1000))
        reads_retrieved += window_reads
        if position[0] > reference_length:
            keep_searching = False
        print(coverages)
        print(len(reads_retrieved))
    return reads_retrieved, coverages


            



