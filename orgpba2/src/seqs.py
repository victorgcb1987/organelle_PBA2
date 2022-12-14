import gzip

from collections import Counter
from distutils.command.config import config
from random import shuffle
from subprocess import run
from tempfile import NamedTemporaryFile

from Bio import SeqIO

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import COMPLEMENTARY_NUCLEOTIDE as compl_nucl
from src.config import OUTPUT_FOLDERS as out_dir
from src.dependencies import get_executables
from src.utils import file_is_compressed, folder_exists, file_exists, sequence_kind


def count_seqs(path):
    #count number of seqs in file
    #by counting number of headers
    if file_is_compressed(path):
        grep_command = ["zgrep"]
    else:
        grep_command = ["grep"]

    if sequence_kind(path) == "fastq":
        grep_command += ["\"@\"", str(path), 
                         "|", "wc -l"]
    elif sequence_kind(path) == "fasta":
        grep_command += ["\">\"", str(path),
                         "|", "wc -l"]
    run_grep_command = run("\t".join(grep_command),
                           shell=True, capture_output=True)
    return int(run_grep_command.stdout.decode())


def get_seqs_id_from_paf_file(paf_file, min_identity=0.7):
    # Recover read names mapped against reference genome
    # Filters by identity (default 70%)
    list_of_seqs = []
    with open(paf_file) as paf_fhand:
        for line in paf_fhand:
            line = line.split()
            seq_name = line[0]
            for item in line:
                if item.startswith("de:f:"):
                    alignment_indentity = (float(1) - float(item.split(":")[2]))
                    if alignment_indentity >= min_identity:
                        if seq_name not in list_of_seqs:
                            list_of_seqs.append(seq_name)
    cmd = ["cut -f1", "{}".format(str(paf_file)), "|", "sort", "|", "uniq"]
    filter_run = run(" ".join(cmd), shell=True, capture_output=True)
    list_of_seqs = filter_run.stdout.decode().rstrip().split("\n")
    return list_of_seqs


def get_reads_total_nucleotides(path):
    if file_is_compressed(path):
        grep_command = "zcat"
    else:
        grep_command = "cat"
    cmd = [grep_command, str(path), "|egrep -e '([ATGCN])\w+'|tr -d '\n'|wc -c"]
    nucleotide_count_run = run(" ".join(cmd), shell=True, capture_output=True)
    gigabase = 1000000000

    return float((nucleotide_count_run.stdout.decode())) / gigabase
    

def store_reads_by_read_name(reads):
    #Sometimes reads have more than one match
    #This function groups matches from the same read
    grouped_reads_matches = {}
    for read in open(reads):
        read = read.split()
        read_name = read[0]
        read_length = int(read[1])
        read_start = int(read[2])
        read_end = int(read[3])
        aln_start = int(read[7])
        aln_end = int(read[8])
        strand = read[4]
        num_matched_nucls = int(read[9])
        aln_length = int(read[10])
        if read_name not in grouped_reads_matches:
            grouped_reads_matches[read_name] = {"read_length": read_length, "hits": [{"read_start": read_start, 
                                                                     "read_end": read_end,
                                                                     "aln_start": aln_start,
                                                                     "aln_end": aln_end,
                                                                     "num_matched_nucls": num_matched_nucls,
                                                                     "strand": strand,
                                                                     "read_length": read_length,
                                                                     "aln_length": aln_length}],
                                                "status": "to check"}
        else:
            grouped_reads_matches[read_name]["hits"].append({"read_start": read_start, 
                                             "read_end": read_end,
                                             "aln_start": aln_start,
                                             "aln_end": aln_end,
                                             "num_matched_nucls": num_matched_nucls,
                                             "strand": strand,
                                             "read_length": read_length,
                                             "aln_length": aln_length})
    return grouped_reads_matches

def mark_reads_as_invalid(grouped_reads_matches, distance_from_assembly_start=1000, matches_cutoff=0.5):
    #criteria
    #it's a gapped alignment
    #different regions of the read maps to overlapping positions
    for read_name, items in grouped_reads_matches.items():
        hits = items["hits"]
        if len(hits) == 1:
            matched_nucls_fraction = float(hits[0]["aln_length"] / hits[0]["read_length"])
            if matched_nucls_fraction >= matches_cutoff:
                grouped_reads_matches[read_name]["status"] = "OK"
            else:
                grouped_reads_matches[read_name]["status"] = "invalid"
            grouped_reads_matches[read_name]["status"] = "OK"
        if len(hits) == 2:
            aln_status = "invalid"
            strands = []
            hit_from_start = False
            for hit in hits:
                strands.append(hit["strand"])
                if hit["aln_start"] <= distance_from_assembly_start:
                    hit_from_start = True
            if hit_from_start and strands[0] == strands[1]:
                aln_status = "OK"
            grouped_reads_matches[read_name]["status"] = aln_status
        if len(hits) > 2:
            grouped_reads_matches[read_name]["status"] = "invalid"
    return grouped_reads_matches

def get_filtered_seqs_id_from_paf_file(paf_file):
    #Return a list of mapped reads that aren't marked as invalid (see mark_reads_as_invalid)
    valid_reads_ids = []
    grouped_read_matches = store_reads_by_read_name(paf_file)
    classified_reads = mark_reads_as_invalid(grouped_read_matches)
    for read_name, info in classified_reads.items():
        if info["status"] == "OK":
            valid_reads_ids.append(read_name)
    return valid_reads_ids, classified_reads

def write_seqs_from_seqs_id(seq_ids, seqs_pool, seqs_out_fpath, overwrite=False):
    #Extract reads mapped to reference genome
    seqtk_exectuable = get_executables(exec_reqs["seqtk"])
    seqs_ids_fhand = NamedTemporaryFile()
    seqs_ids_fhand.write("\n".join(seq_ids).encode())
    seqs_ids_fhand.flush()
    if not file_exists(seqs_out_fpath) or overwrite:
        cmd = [seqtk_exectuable, "subseq", str(seqs_pool), 
               seqs_ids_fhand.name, ">", str(seqs_out_fpath)]
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

def get_seq_length(seq_fpath):
    #Get number of nucleotides from single record fasta
    if file_is_compressed(seq_fpath):
        seq_fhand = gzip.open(seq_fpath, "rt")
    else:
        seq_fhand = open(seq_fpath, "rt")
    record = SeqIO.read(seq_fhand, "fasta")
    nucleotide_counts = Counter(record.seq)
    num_nucleotides = 0
    for nucl, number in nucleotide_counts.items():
        num_nucleotides += number
    seq_fhand.close()
    return num_nucleotides

def get_seqs_length(seqs_fpath):
    #Get number of nucleotides from multi record fasta
    records = SeqIO.parse(seqs_fpath, "fasta")
    num_nucleotides = 0
    for record in records:
        nucleotide_counts = Counter(record.seq)
        for nucl, number in nucleotide_counts.items():
            num_nucleotides += number
    return num_nucleotides

def reverse_complement(sequence):
    #Create reverse complement sequence from a original sequence
    return "".join([compl_nucl[nucleotide.upper()] for nucleotide in sequence[::-1]])

def complement(sequence):
    #Create complement sequence from a original sequence
    return "".join([compl_nucl[nucleotide.upper()] for nucleotide in sequence])


def find_origin(blastn_input, margin=10):
    #Find origin from an assembly by blast between assembly and reference genome
    #Assumes that distance from origin is equal or less than the margin
    origin_found = False
    #blastn_input must be formatted with -outfmt '6 std qlen slen'
    with open(blastn_input) as blastn:
        for line in blastn:
            line = line.split()
            reference_start = int(line[8])
            reference_end = int(line[9])
            query_start = int(line[6])
            query_end = int(line[7])
            if reference_start >= 1 and reference_start <= 1+margin:
                origin_found = {"reference": reference_start-1, 
                                "query": query_start-1, 
                                "strand": "+"}
            elif reference_end >= 1 and reference_end <= 1+margin:
                origin_found = {"reference": reference_end-1, 
                                "query": query_end, 
                                "strand": "-"}
            if origin_found:
                return origin_found
        return origin_found


def run_blast(query, database, options, blast_arguments=None):
    # Runs blast and returns results
    if blast_arguments is None:
        msg += "No blast options found"
        raise ValueError(msg)
    else:
        out_fpath = options["out_dir"] / out_dir["blast_db"]
        if not folder_exists(out_fpath):
            out_fpath.mkdir(parents=True, exist_ok=True)
        makeblastdb =  get_executables(exec_reqs["makeblastdb"])
        if blast_arguments["kind"] == "blastn" or blast_arguments["kind"] == "tblastn":
            kind = "nucl"
        else:
            kind = "prot"
        dbtype = "-dbtype {}".format(kind)
        dbname = "-in {}".format(str(database))
        cmd = [str(makeblastdb), dbtype, dbname]
        makeblastdb_run = run(" ".join(cmd), shell=True, 
                              capture_output=True)
        if makeblastdb_run.returncode != 0:
            raise RuntimeError(makeblastdb_run.stderr.decode())
        cmd = [blast_arguments["kind"]]
        if "task" in blast_arguments:
            cmd.append("-task {}".format(blast_arguments["task"]))
        query_cmd = "-query {}".format(str(query))
        database_cmd = "-db {}".format(str(database))
        cmd.append(query_cmd)
        cmd.append(database_cmd)
        if "outmft" in blast_arguments:
            cmd.append(blast_arguments["outmft"])
        output_fname = "{}_against_{}.{}".format(query.stem, 
                                                 database.stem, 
                                                 blast_arguments["kind"])
        blast_output_fpath = out_fpath / output_fname
        cmd.append("> {}".format(str(blast_output_fpath)))
        blast_run = run(" ".join(cmd), shell=True, 
                        capture_output=True)
        results = {"output_file": blast_output_fpath,
                   "return_code": blast_run.returncode,
                   "log_messages": blast_run.stderr.decode()}
        return results

def find_circularity(assembly, options, margin_length=10, overlap_length=100):
    #Finds assembly's circularty by self blast
    query_length = get_seq_length(assembly)
    blast_arguments = {"kind": "blastn",
                        "task": "megablast",
                        "outmft": "-outfmt '6 std qlen slen'" }
    blastn_results = run_blast(assembly, assembly, options, blast_arguments=blast_arguments)
    blast_results_output = blastn_results["output_file"]
    overlaps = False
    with open(blast_results_output) as blast_ouput:
        found_overlap_start = False
        found_overlap_end = False
        for line in blast_ouput:
            line = line.rstrip().split()
            reference_start = int(line[8])
            reference_end = int(line[9])
            query_start = int(line[6])
            query_end = int(line[7])
            hit_identity = float(line[2])
            if query_start == 1 and query_end == query_length and reference_start == 1 and reference_end == query_length:
                continue
            else:
                if query_start <= margin_length and query_end >= (overlap_length + query_start):
                    if not found_overlap_start:
                        found_overlap_start = True
                        begin_overlap_start = query_start
                        begin_overlap_end = query_end
                elif query_start <= (query_length - overlap_length) and query_end >= (query_length - margin_length):
                    if not found_overlap_end:
                        found_overlap_end = True
                        end_overlap_start = query_start
                        end_overlap_end = query_end
                if reference_start <= margin_length and reference_end >= (reference_start + overlap_length):
                    if found_overlap_start:
                        found_overlap_start = True
                        begin_overlap_start = reference_start
                        begin_overlap_end = reference_end
                elif reference_start <= (query_length - overlap_length) and reference_end >= (query_length - margin_length):
                    if found_overlap_end:
                        found_overlap_end = True
                        end_overlap_start = reference_start
                        end_overlap_end = reference_end
            if found_overlap_start and found_overlap_end:
                overlaps = {"overlap_start": (begin_overlap_start -1, begin_overlap_end -1),
                            "overlap_end": (end_overlap_start -1 , end_overlap_end -1 ),
                            "overlap_identity": hit_identity}
        
    return overlaps


def remove_circularity_redundancy(assembly, circularity):
    #removes circularity found in find_circularity function
    assembly_fhand = open(assembly)
    header = assembly_fhand.readline()
    no_redundant_seq = header
    seq = "".join([line.rstrip() for line in assembly_fhand])
    seq = seq[:circularity["overlap_end"][0]]
    no_redundant_seq += seq
    assembly_fhand.close()
    return no_redundant_seq


def reconstruct_assembly_from_origin(sequence, origin):
    #Reconstructs assembly if reference genome origin was found in assembly
    sequence = "".join(sequence.split("\n")[1:])
    if origin["strand"] == "-":
        end_sequence = sequence[:origin["query"]]
        origin_sequence = sequence[origin["query"]:]
        reconstructed_sequence = origin_sequence+end_sequence
        reconstructed_sequence = reverse_complement(reconstructed_sequence)
    elif origin["strand"] =="+":
        origin_sequence = sequence[origin["query"]:]
        end_sequence = sequence[:origin["query"]]
        reconstructed_sequence = origin_sequence+end_sequence
        
    return reconstructed_sequence


def calculate_repeats_regions(nucmer_output):
    #Find IRs regions in assembly by using nucmer results between reference and assembly
    with open(nucmer_output) as coordinates_fhand:
        begin_matches = False
        matches = {}
        for line in coordinates_fhand:
            if line.startswith("="):
                begin_matches = True
                continue
            if not begin_matches:
                continue
            ref_start, ref_end, _, query_start, query_end, _, ref_length, query_length, _, identity, _, ref_id, assembly_id = line.split()
            
            ref_start = int(ref_start)
            ref_end = int(ref_end)
            query_start = int(query_start)
            query_end = int(query_end)
            query_length = int(query_length)
            
            if ref_start > ref_end:
                ref_start, ref_end = ref_end, ref_start
            else:
                ref_start, ref_end = ref_start, ref_end
            if query_start > query_end:
                query_start, query_end = query_end, query_start
            else:
                query_start, query_end = query_start, query_end
            if query_length in matches:
                repeat = True
                regions = [(ref_start, ref_end), (query_start, query_end)]
                for region in regions:
                    if region not in matches[query_length]:
                        repeat = False
                if repeat:
                    return sorted(regions)
            else:
                matches[query_length] = [(ref_start, ref_end), (query_start, query_end)]
    return None


def find_biggest_inverted_repeat_sequence(inverted_repeat_regions, nucmer_output):
    #Selects the largest IR region found in assembly
    with open(nucmer_output) as coordinates_fhand:
        begin_matches = False
        best_match = {}
        for line in coordinates_fhand:
            if line.startswith("="):
                begin_matches = True
                continue
            if not begin_matches:
                continue
            ref_start, ref_end, _, query_start, query_end, _, ref_length, query_length, _, identity, _, ref_id, query_id = line.split()
            ref_start = int(ref_start)
            ref_end = int(ref_end)
            query_start = int(query_start)
            query_end = int(query_end)
            query_length = int(query_length)
            for inverted_repeat in inverted_repeat_regions:
                repeat_start = inverted_repeat[0]
                repeat_end = inverted_repeat[1]
                if ref_start > repeat_end or ref_end < repeat_start:
                    continue
                if ref_start < repeat_start:
                    ref_start = repeat_start
                if ref_end > repeat_end:
                    ref_end = repeat_end
                if query_start > query_end:
                    query_start, query_end = query_end, query_start
                if not best_match:
                    best_match = {"contig_name": query_id, "length": ref_end - ref_start, "contig_position": [query_start, query_end],
                                  "ref_position": [ref_start, ref_end]}
                elif best_match["contig_name"] == query_id:
                    if ref_start < best_match["ref_position"][-1]:
                        continue
                    distance = abs(best_match["contig_position"][-1] - query_start)
                    if distance < 500:
                        best_match["length"] += ref_end - ref_start
                        best_match["contig_position"] = [best_match["contig_position"][0], query_end]
                        best_match["ref_position"] = [best_match["ref_position"][0], ref_end]
                else:
                    if (ref_end - ref_start) > best_match["length"]:
                        best_match = {"contig_name": query_id, "length": ref_end - ref_start, "contig_position": [query_start, query_end], 
                                      "ref_position": [ref_start, ref_end]}
                        break
    return best_match


def check_if_assembly_is_complete(assembly_fpath, genome_size, margin=10):
    #Checks if assembly is complete
    #If more than one contig is in the fasta file, it will select the largest one for next steps
    needs_scaffolding_step = True
    number_of_contigs = count_seqs(assembly_fpath)
    if number_of_contigs == 1:
        print("only one contig generated, no scaffolding")
        record = SeqIO.read(assembly_fpath, "fasta")
        if len(record.seq) >= genome_size - margin:
            needs_scaffolding_step = False
        return needs_scaffolding_step, record
    records = SeqIO.parse(assembly_fpath, "fasta")
    max_length = {"id": "record", "length": 0}
    for record in records:
        if len(record.seq) > max_length["length"]:
            max_length["record"] = record
            max_length["length"] = len(record.seq)
        if len(record.seq) >= genome_size - margin:
            needs_scaffolding_step = False
    return needs_scaffolding_step, max_length["record"]

def create_renamed_assembly(assembly_fpath, out_fhand):
    #Renames contigs by contig length
    records = SeqIO.parse(assembly_fpath, "fasta")
    num_contig = 1
    for record in records:
        length = len(record.seq)
        out_fhand.write(">contig{}_size{}\n".format(str(num_contig), str(length)))
        out_fhand.write(str(record.seq)+"\n")
        out_fhand.flush()
        num_contig += 1
    out_fhand.close()


def filter_overlapping_contigs(assembly, options, filtered_fpath, overlap_fraction=0.8):
    #Sometimes contigs can share a large margin of the reference assmebly
    #It filters away overlapping contigs, keeping the largest one
    #Do a selfblastn
    #check if contigs have at least a 0.8 of overlapping
    #then keep the largest
    blast_arguments = {"kind": "blastn",
                        "task": "megablast",
                        "outmft": "-outfmt '6 std qlen slen'" }
    blastn_results = run_blast(assembly, assembly, options, blast_arguments=blast_arguments)

    overlapped_seqs = []
    with open(blastn_results["output_file"]) as blast_ouput:
        for line in blast_ouput:
            line = line.split()
            query_name = line[0]
            subject_name = line[1]
            alignment_length = int(line[3])
            query_length = int(line[-2])
            if query_name == subject_name:
                continue
            else:
                if float(alignment_length / query_length ) >= overlap_fraction and query_name not in overlapped_seqs:
                    overlapped_seqs.append(query_name)
    records = SeqIO.parse(assembly, "fasta")
    out_fhand = open(filtered_fpath, "w")
    for record in records:
        if record.id in overlapped_seqs:
            continue
        else:
            SeqIO.write(record, out_fhand, "fasta")
    out_fhand.flush()
    out_fhand.close()


def write_circular_region(redundant_seq_fpath, circularity, output_fpath):
    circular_region_start_0 = circularity["overlap_start"][0]
    circular_region_start_1 = circularity["overlap_start"][1]
    circular_region_end_0 = circularity["overlap_end"][0]
    circular_region_end_1 = circularity["overlap_end"][1]
    print(circular_region_start_0, circular_region_start_1)
    header = ">circular_sequence, found at {}-{} and {}-{}\n".format(circular_region_start_0, circular_region_start_1,
                                                                    circular_region_end_0, circular_region_end_1)
    record = SeqIO.read(redundant_seq_fpath, "fasta")
    circular_sequence = record.seq[circular_region_start_0: circular_region_start_1]
    with open(output_fpath, "w") as out_fhand:
        out_fhand.write(header)
        out_fhand.write(str(circular_sequence))
        out_fhand.flush()
    

    




    
        