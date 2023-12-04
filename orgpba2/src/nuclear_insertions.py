from statistics import median

from Bio import SeqIO

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FOLDERS as out_dir
from src.dependencies import get_executables
from tempfile import NamedTemporaryFile
from src.utils import file_exists
from subprocess import run


def get_reads_alignments_info(reads_fhand, organelle_length=0, repeats=False, exclude_potential_chimeras=True):
#Get info of diferent reads from an archive in paf format
    if repeats:
        offset = repeats[1][0] - repeats[0][0]
    else:
        offset = False
    print(offset)
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
            if subject_start >= organelle_length and subject_end >= organelle_length:
                subject_start = subject_start - organelle_length
                subject_end = subject_end - organelle_length
            if offset:
                print(subject_start, subject_end)
                if subject_start >= repeats[1][0] and subject_end <= repeats[1][1]:
                    subject_start = (subject_start - repeats[1][0]) + repeats[0][0]
                    subject_end = (repeats[1][1] - subject_end) - repeats[0][1] 
                    print("YES", subject_start, subject_end)
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
            elif strand != reads_alignments_info[read_name]["strand"] and exclude_potential_chimeras:
                reads_alignments_info.pop(read_name)
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


def write_seqs_from_seqs_id(seq_ids, seqs_pool, seqs_out_fpath, overwrite=False):
    #Extract reads mapped to reference genome
    seqtk_exectuable = get_executables(exec_reqs["seqtk"])
    seqs_ids_fhand = NamedTemporaryFile()
    seqs_ids_fhand.write("\n".join(seq_ids).encode())
    seqs_ids_fhand.flush()
    if not file_exists(seqs_out_fpath) or overwrite:
        cmd = [seqtk_exectuable, "subseq", str(seqs_pool), 
               seqs_ids_fhand.name, "| gzip -c >", str(seqs_out_fpath)]
        seqtk_run  = run(" ".join(cmd), shell=True,
                         capture_output=True)
        results = {"output_files": [seqs_out_fpath],
                   "return_code": seqtk_run.returncode,
                   "log_messages": seqtk_run.stderr.decode()}
    else:
        results = {"output_files": [seqs_out_fpath],
                   "return_code": 0,
                   "log_messages": "File exists"}
    return results


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


def remove_organelle_offset(alignments_info, organelle_length):
    for readname, value in alignments_info.items():
        if value['subject_start'] >= organelle_length and value['subject_end'] >= organelle_length:
            alignments_info[readname]['subject_start'] = value['subject_start'] - organelle_length
            alignments_info[readname]['subject_end'] = value['subject_end'] - organelle_length
    return alignments_info


def exclude_reads_with_less_coverage(organelles_alignments_info, exclude_alignments_info):
    alignments_info = {}
    for read_name, values in organelles_alignments_info.items():
        if not exclude_alignments_info.get(read_name, False):
            alignments_info[read_name] = values
        elif values['total_aligned_length'] >= exclude_alignments_info[read_name]['total_aligned_length']:
            alignments_info[read_name] = values
    return alignments_info


def group_reads_of_same_insertion(insertion_reads, organelle_boundaires=1000):
    # insertions_positions[name] = {'insertion_starts' : joined_groups_starts, 
    #                               'insertion_ends' : joined_groups_end,
    #                                'organelle_start' : p_start_organ, 
    #                                'organelle_end' : p_end_organ,
    #                                'nuclear' : chrom,
    #                                'reads' : list_reads}
    groups = []
    # Recorrer el diccionario original y agrupar los valores de 'name' según si las regiones se solapan o no
    for key, value in insertion_reads.items():
        group_found = False
    
    # Buscar si existe algún grupo que se solape con esta entrada
        for group in groups:
            if not groups:
                break
            group_start = min(group['insertion_starts'])
            group_end = max(group['insertion_ends'])
            group_organelle_start = min(group['organelle_starts'])
            group_organelle_end = max(group['organelle_ends'])
            if value["chrom"] == group["nuclear"]:
                if abs(value['organelle_start'] - group_organelle_start) <= organelle_boundaires or abs(value['organelle_end'] - group_organelle_end) <= organelle_boundaires:
                    if abs(value['insertion_start'] - group_start) <= organelle_boundaires or abs(value['insertion_end'] - group_end) <= organelle_boundaires:
                        group['readnames'].append(key)
                        group["insertion_starts"].append(value['insertion_start'])
                        group["insertion_ends"].append(value["insertion_end"])
                        group["organelle_starts"].append(value["organelle_start"])
                        group["organelle_ends"].append(value["organelle_end"])
                        group_found = True
                        break

    # Si no se encontró ningún grupo que se solape, crear uno nuevo
        if not group_found:
            groups.append({
                'organelle_starts': [value['organelle_start']],
                'organelle_ends': [value['organelle_end']],
                'readnames': [key],
                'nuclear': value["chrom"],
                'insertion_starts': [value["insertion_start"]],
                'insertion_ends': [value["insertion_end"]]
            })
    return groups
