import argparse
from pathlib import Path

from src.curation import run_nucmer, get_genomic_coordinates_from_assembly, create_curated_genome


from src.seqs import get_seq_length
from src.utils import check_results


def parse_arguments():
    desc = "Curate chloroplast assembly by homology with a reference genome"
    parser = argparse.ArgumentParser(description=desc)
    help_reference = "(Required) Reference genome."
    parser.add_argument("--reference_sequence", 
                        "-r", type=str,
                        help=help_reference,
                        required=True)
    help_reference = "(Required) Assembly to curate."
    parser.add_argument("--assembly", 
                        "-a", type=str,
                        help=help_reference,
                        required=True)
    help_output = "(Required) Output dir"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    return parser

#Parse and return values given to options when running this program
def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    reference_fpath = Path(options.reference_sequence)
    assembly_fpath = Path(options.assembly)
    output_dir = Path(options.out)
    
    return {'reference_input': reference_fpath,
            'assembly_fpath': assembly_fpath, 
            'out_dir': output_dir}

def main():
    options = get_options()
    msg = "Curating genome assenbly with nucmer"
    nucmer_results = run_nucmer(options, overwrite=True)
    check_results(msg, nucmer_results)
    nucmer_output = nucmer_results["output_file"]
    genomic_coordinates = get_genomic_coordinates_from_assembly(nucmer_output)
    curated_assembly_results = create_curated_genome(options, genomic_coordinates)
    curated_assembly = curated_assembly_results["output_dir"] / "curated_assembly.fasta"
    curated_assembly_length = get_seq_length(curated_assembly)
    print("Length of curated assembly: {}\n".format(curated_assembly_length))


if __name__ == "__main__":
    main()