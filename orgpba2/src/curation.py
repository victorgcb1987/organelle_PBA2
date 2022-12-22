import shutil

from itertools import groupby
from operator import itemgetter
from pathlib import Path
from subprocess import run

from Bio import SeqIO

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.config import OUTPUT_FOLDERS as out_dir
from src.config import OUTPUT_FILENAMES as outfile
from src.dependencies import get_executables
from src.utils import folder_exists, chunkstring

def run_nucmer(options, overwrite=False):
    assembly_fname = options["assembly_fpath"].stem
    reference_fname = options["reference_input"].stem
    prefix = "nucmer_{}_{}".format(reference_fname, assembly_fname)
    run_nucmer_dir = options["out_dir"] / out_dir["curated_assembly"]
    #check if output folder exists
    if folder_exists(run_nucmer_dir) and overwrite:
        if overwrite:
            shutil.rmtree(run_nucmer_dir)
            run_nucmer_dir.mkdir(parents=True, exist_ok=True) 
        else:
            msg = "output dir {} already exists, skipping"
            msg += " nucmer_alignment"
            results = {"output_files": run_nucmer_dir,
                       "return_code": 0,
                       "log_messages": msg.format(str(run_nucmer_dir))}
            return results
    else:
         run_nucmer_dir.mkdir(parents=True, exist_ok=True) 
    nucmer_executable = get_executables(exec_reqs["nucmer"])
    cmd = [nucmer_executable, "--maxmatch", "--nosimplify", str(options["reference_input"]), str(options["assembly_fpath"]),
          "-p", prefix, "--coords"]
    nucmer_run = run(" ".join(cmd), shell=True, 
                       capture_output=True)

    current_path = Path.cwd()
    for filename in ["{}.coords", "{}.delta"]:
         current_filepath = current_path /filename.format(prefix)
         destiny_filepath = run_nucmer_dir / filename.format(prefix)
         shutil.move(current_filepath, destiny_filepath)
  
    results = {"output_file": run_nucmer_dir / "{}.coords".format(prefix),
               "return_code": nucmer_run.returncode,
               "log_messages": nucmer_run.stderr.decode()}

    return results
    

def get_genomic_coordinates_from_assembly(nucmer_output):
    genomic_coordinates = []
    with open(nucmer_output) as coordinates_fhand:
        begin_matches = False
        for line in coordinates_fhand:
            if line.startswith("="):
                begin_matches = True
                continue
            if not begin_matches:
                continue
            ref_start, ref_end, _, assembly_start, assembly_end, _, ref_length, assembly_length, _, identity, _, ref_id, assembly_id = line.split()
            
            ref_start = int(ref_start)
            ref_end = int(ref_end)
            assembly_start = int(assembly_start)
            assembly_end = int(assembly_end)
            
            if assembly_start > assembly_end:
                assembly_strand = "-"
                assembly_start, assembly_end = assembly_end, assembly_start
            else:
                assembly_strand = "+"
            if assembly_start < 50:
                assembly_start = 1
            block = {"assembly_start": int(assembly_start), "assembly_end": int(assembly_end),
                                 "ref_start": int(ref_start), "ref_end": int(ref_end),
                                  "contig_id": assembly_id, "strand": assembly_strand} 
            if not genomic_coordinates:
                genomic_coordinates.append(block)
                continue
            last_block = genomic_coordinates[-1]

            if ref_start >= last_block["ref_end"]:
                genomic_coordinates.append(block)
                
            elif ref_start < last_block["ref_end"] and ref_start >= last_block["ref_start"] and ref_end >= last_block["ref_end"]:
                offset = ref_end - last_block["ref_end"]
                if block["strand"] == "+":
                   block["assembly_start"] = block["assembly_end"] - offset + 1
                if block["strand"] == "-":
                   block["assembly_end"] = block["assembly_start"] + offset - 1
                genomic_coordinates.append(block)
    return genomic_coordinates


def create_curated_genome(options, coordinates_blocks):
    assembly_fpath = str(options["assembly_fpath"])
    curated_assembly_fpath = options["out_dir"] / out_dir["curated_assembly"] / outfile["curated_assembly"]
    curated_contigs_fpath = options["out_dir"] / out_dir["curated_assembly"] / outfile["curated_contigs"]
    assembly_records = SeqIO.index(assembly_fpath, "fasta")
    output_fhand = open(curated_assembly_fpath, "w")
    output_contigs_fhand = open(curated_contigs_fpath, "w")
    seq = ""
    contig_count = 0
    for block in coordinates_blocks:
        if abs(block["assembly_start"] - block["assembly_end"]) < 10:
            continue
        contig_count += 1
        if block["strand"] == "+":
            subseq = str(assembly_records[block["contig_id"]].seq[block["assembly_start"] -1 : block["assembly_end"]])
            output_contigs_fhand.write(">curated_contig{}\t{}\n{}\n".format(contig_count, len(subseq), str(subseq)))
            seq += subseq
        if block["strand"] == "-":
            subseq = str(assembly_records[block["contig_id"]].seq[block["assembly_start"] -1 : block["assembly_end"]].reverse_complement())
            output_contigs_fhand.write(">curated_contig{}\t{}\n{}\n".format(contig_count, len(subseq), str(subseq)))
            seq += subseq
            total_length = len(seq)
    seq = chunkstring(seq, 60)
    output_fhand.write(">curated_assembly\t{}\n".format(total_length))
    for line in seq:
        output_fhand.write(line)
        output_fhand.flush()
    results = {"output_dir": options["out_dir"] / out_dir["curated_assembly"]}
    return results


def retrieve_not_overlapping_coordinates_from_contigs(contig_with_largest_ir, nucmer_output):
    coordinates = {}
    matches = {}
    with open(nucmer_output) as coordinates_fhand:
        begin_matches = False
        for line in coordinates_fhand:
            if line.startswith("="):
                begin_matches = True
                continue
            if not begin_matches:
                continue
            ref_start, ref_end, _, query_start, query_end, _, ref_length, assembly_length, _, identity, _, ref_id, query_id = line.split()
            ref_start = int(ref_start)
            ref_end = int(ref_end)
            query_start = int(query_start)
            query_end = int(query_end)
            if query_start > query_end:
                query_start, query_end = query_end, query_start
            if ref_start == query_start and ref_end == query_end and ref_id == query_id:
                coordinates[query_id] = range(query_start, query_end)
                if query_id == contig_with_largest_ir["contig_name"]:
                    ir_start, ir_end = query_start, query_end
            else: 
                if ref_id not in matches:
                    matches[ref_id] = [range(ref_start, ref_end)]
                else:
                    matches[ref_id].append(range(ref_start, ref_end))
                if query_id not in matches:
                    matches[query_id] = [range(query_start, query_end)]
                else:
                    matches[query_id].append(range(query_start, query_end))
    
    regions_without_overlapping = {}
    for contig, coordinate in coordinates.items():
        ranges = []
        coordinate_set = set(coordinate)
        if contig not in matches:
            regions_without_overlapping[contig] = [(coordinate[0], coordinate[-1])]
        else:
            for overlap in matches[contig]:
                coordinate_set = coordinate_set.difference(overlap)
        
            for k, g in groupby(enumerate(coordinate_set), lambda x:x[0]-x[1]):
                group = (map(itemgetter(1), g))
                group = list(map(int,group))
                ranges.append((group[0]+1, group[-1]))
            regions_without_overlapping[contig] = ranges
    print(regions_without_overlapping)
    print(contig_with_largest_ir['contig_position'][0])
     #Aqui hay varias posibilidades:
    #1 Los contigs no tienen solapamiento, de modo que tienen el mismo tamaño
    #2 No hay solapamiento, pero la detección del ir se ha hecho con una especie distante y por lo tanto la longitud detectada es menor
    #3 Los contigs tiene solapmiento y al final desparece el ir porque todo solapa
    if not regions_without_overlapping[contig_with_largest_ir["contig_name"]]:
        regions_without_overlapping[contig_with_largest_ir["contig_name"]] = [(ir_start, ir_end)]
    else:
        not_overlapping_contig_length = abs(regions_without_overlapping[contig_with_largest_ir["contig_name"]][0][0] - regions_without_overlapping[contig_with_largest_ir["contig_name"]][0][1])
        ir_detected_length = abs(contig_with_largest_ir['contig_position'][0] - contig_with_largest_ir['contig_position'][1])
        if not_overlapping_contig_length < ir_detected_length :
            regions_without_overlapping[contig_with_largest_ir["contig_name"]] = [(contig_with_largest_ir['contig_position'][0], contig_with_largest_ir['contig_position'][-1])]
    print("Not overlapping", regions_without_overlapping)
    return regions_without_overlapping

def write_not_overlapping_regions(options, not_overlapping_regions):
    asssembled_contigs_fpath = options["out_dir"] / out_dir["flye"] / outfile["assembly"]
    contigs_without_overlaps = options["out_dir"] / out_dir["curated_assembly"] / outfile["contigs_without_overlaps"]
    print(contigs_without_overlaps)
    assembly_input = SeqIO.index(str(asssembled_contigs_fpath), "fasta")
    with open(str(contigs_without_overlaps), "w") as out_fhand:
        for contig, regions in not_overlapping_regions.items():
            for region in regions:
                print(regions)
                if len(regions) == 1:
                    subseq = assembly_input[contig][region[0]-1:region[-1]]
                    subseq.id = "{}_{}-{}".format(contig, region[0], region[-1])
                else:
                    subseq = assembly_input[contig][region[0]-1:region[-1]]
                    subseq.id = "{}_{}-{}".format(contig, region[0], region[-1])
                SeqIO.write(subseq, out_fhand, "fasta")
                out_fhand.flush()
    return contigs_without_overlaps


    


def orient_and_concatenate_contigs(options, nucmer_output, min_identity=90):
    contigs_fpath = options["not_overlapping_contigs"]
    concatenated_out_fpath = options["out_dir"] / out_dir["curated_assembly"] / outfile["concatenated_contigs"]
    oriented_out_fpath = options["out_dir"] / out_dir["curated_assembly"] / outfile["oriented_contigs"]
    matches = []
    with open(nucmer_output) as coordinates_fhand:
        begin_matches = False
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
            identity = float(identity)
            if identity < min_identity:
                continue
            if query_start > query_end:
                query_start, query_end = query_end, query_start
                strand = "-"
            else:
                strand = "+"
            match = {"position": [ref_start, ref_end], "contig": query_id, 
                     "query_length": query_length, "strand": strand, "identity": identity}

            if not matches:
                matches.append(match)
            else:
                previous_match_range_of_coordinates = set(range(matches[-1]["position"][0], matches[-1]["position"][-1]))
                match_range_of_coordinates = set(range(ref_start, ref_end))
                if previous_match_range_of_coordinates & match_range_of_coordinates:
                    if identity > matches[-1]["identity"] and query_length > matches[-1]["query_length"]:
                        matches.pop()
                        matches.append(match)
                elif identity > matches[-1]["identity"] and query_id == matches[-1]["contig"]:
                    matches.pop()
                    matches.append(match)
                elif query_id != matches[-1]["contig"]:
                    matches.append(match)

    print("XXXXXXX", matches)

    contigs = SeqIO.index(str(contigs_fpath), "fasta")
    with open(str(concatenated_out_fpath), "w") as out_fhand:
        seq = ""
        for match in matches:
            if match["strand"] == "+":
                subseq = str(contigs[match["contig"]].seq[:])
            elif match["strand"] == "-":
                subseq = str(contigs[match["contig"]].seq[:].reverse_complement())
            seq += subseq
        out_fhand.write(">Concatenated_contigs\t{}\n{}".format(len(seq), str(seq)))
        out_fhand.flush()
    
    with open(str(oriented_out_fpath), "w") as out_fhand:
        for match in matches:
            if match["strand"] == "+":
                seq = str(contigs[match["contig"]].seq[:])
            elif match["strand"] == "-":
                seq = str(contigs[match["contig"]].seq[:].reverse_complement())
            out_fhand.write(">{} strand:{} {}\n{}\n".format(match["contig"], match["strand"], len(seq), seq))
            out_fhand.flush()


        



