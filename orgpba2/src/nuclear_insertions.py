from Bio import SeqIO

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import COMPLEMENTARY_NUCLEOTIDE as compl_nucl
from src.config import OUTPUT_FOLDERS as out_dir
from src.utils import sequence_kind


def get_reads_alignments_info(reads_fhand):
#Get info of diferent reads from an archive in paf format
    reads_alignments_info = {}
    for line in reads_fhand:
        if line:
            line = line.rstrip()
            line = line.split()
            read_name = line[0]
            read_length = int(line[1])
            strand = line[4]
            subject_name = line[5]
            target_positions = int(line[7]), int(line[8])
            subject_start = int(line[7])
            subject_end = int(line[8])
            query_start = int(line[2])
            query_end = int(line[3])
            strand = line[4]
            total_alignment =  int(line[3]) - int(line[2])
            if read_name not in reads_alignments_info:
                reads_alignments_info[read_name] = {'length' : read_length, 
                                                    'target_positions' : target_positions,
                                                    'total_aligned_length': total_alignment,
                                                    'subject_start': subject_start,
                                                    'subject_end': subject_end,
                                                    'query_start': query_start,
                                                    'query_end': query_end,
                                                    'strand': strand,
                                                    'subject_name': subject_name}
    return reads_alignments_info

def calculate_reads_query_coverage(alignments_info):
#Calculate reads query coverage (aligned part respect total)
    coverages = {}
    for read, alignment in alignments_info.items():
        read_length = alignment['length']
        coverage = float(alignment['total_aligned_length'] / read_length)
        coverages[read] = coverage
    return coverages

def select_reads_by_coverage(reads_coverage, coverage_cutoff = 0.9, mode="under"):
#Obtain reads that their query coverage is under a determinate acceptation value
    if mode == "under":
        return [read_name for read_name, coverage in reads_coverage.items() if coverage <= coverage_cutoff]
    elif mode == "over":
        return [read_name for read_name, coverage in reads_coverage.items() if coverage >= coverage_cutoff]

def filter_reads_by_name(alignments_info, readnames):
    return {readname: alignments_info[readname] for readname in readnames}

def write_aligned_reads_into_fasta(arguments,filtered_reads,not_aligned_file):
    sequences = arguments["sequences"]
    sequences_kind = sequence_kind(sequences)
    if sequences_kind == "fasta":
        sequences_records = SeqIO.index(str(sequences), "fasta")
    elif sequences_kind == "fastq":
        sequences_records = SeqIO.index(str(sequences), "fastq")
    with open(not_aligned_file, "w") as out_fhand:
        for read in filtered_reads:
            sequence = "{} \n".format(sequences_records[read].seq)
            header = ">{} \n".format(sequences_records[read].id)
            out_fhand.write(header)
            out_fhand.write(sequence)
            out_fhand.flush()
    return not_aligned_file


def get_insertions_positions(nuclear_alignments, organelle_alignments):
        #FEATURES -> COSAS DEL CLOROPLASTOS 
        #NUCLEAR HITS -> COSAS DEL NUCLE
        chroms = []
        insertion_reads = {}
        for readname, features in organelle_alignments.items():
            if readname in nuclear_alignments:
                nuclear_hit = nuclear_alignments[readname]
                insertion_start = 0
                insertion_end = 0
                organelle_start = features["subject_start"]
                organelle_end = features["subject_end"]
                if features["strand"] == "-" and nuclear_hit["strand"] == "-":
                    insertion_start = nuclear_hit["subject_start"] + features["length"] - features["query_end"]
                    insertion_end = insertion_start + (features["query_end"] - features["query_start"])
                    
                    #print(readname, features["strand"], nuclear_hit["strand"], insertion_start, insertion_end)
                if features["strand"] == "+" and nuclear_hit["strand"] == "+":
                    insertion_start = nuclear_hit["subject_start"] + (features["query_start"])
                    insertion_end = insertion_start + ((features["query_end"] - features["query_start"]))
                    
                    #print(readname, features["strand"], nuclear_hit["strand"], insertion_start, insertion_end)
                if features["strand"] == "+" and nuclear_hit["strand"] == "-":
                     insertion_start = nuclear_hit["subject_start"] + (features["length"] - features["query_end"])
                     insertion_end = insertion_start + (features["query_end"] - features["query_start"]) 
                     
                if features["strand"] == "-" and nuclear_hit["strand"] == "+":
                     insertion_end = nuclear_hit["subject_start"] + (features["query_end"])
                     insertion_start = insertion_end - (features["query_end"] - features["query_start"])
                     
                # name_to_insert = readname
                # if name_to_insert in insertion_reads:
                #     name_to_insert = "{}_1".format(readname[:-2])
                insertion_reads[readname]= {'read' : readname, 
                                       'insertion_start' : insertion_start, 
                                       'insertion_end' : insertion_end, 
                                       'organelle_start' : organelle_start, 
                                       'organelle_end': organelle_end,
                                       'chrom': nuclear_hit["subject_name"],
                                       'nuclear_strand': nuclear_hit["strand"],
                                       'organelle_strand': features["strand"]}
                if nuclear_hit["subject_name"] not in chroms:
                    chroms.append(nuclear_hit["subject_name"])
        return chroms, insertion_reads

def exclude_reads_with_less_coverage(organelles_alignments_info, exclude_alignments_info):
    alignments_info = {}
    for read_name, values in organelles_alignments_info.items():
        if not exclude_alignments_info.get(read_name, False):
            alignments_info[read_name] = values
        elif values['total_aligned_length'] >= exclude_alignments_info[read_name]['total_aligned_length']:
            alignments_info[read_name] = values
    return alignments_info
