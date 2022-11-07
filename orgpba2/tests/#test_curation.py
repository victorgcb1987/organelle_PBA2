import unittest

from shutil import rmtree as remove_folder
from pathlib import Path

from src.seqs import calculate_repeats_regions, find_biggest_inverted_repeat_sequence
from src.curation import (run_nucmer, get_genomic_coordinates_from_assembly, create_curated_genome, 
                          retrieve_not_overlapping_coordinates_from_contigs, write_not_overlapping_regions,
                          orient_and_concatenate_contigs)

class TestCuration(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.reference_sequence = test_path / "artha_chloro_ref.fasta"
        self.contigs = test_path / "artha_chloro_parts.fasta"
        self.overlapping_contigs = test_path / "artha_chloro_parts.overlapping.fasta"
        self.run_project_out_dir = test_path / "test_out"
        self.out_dir = self.run_project_out_dir / "04_curated_assembly_mummer"
        self.not_overlapping_regions_seq = self.out_dir / "not_overlapping_seqs.fasta"
        self.oriented_contigs = self.out_dir / "oriented_contigs.fasta"
        self.concatenated_contigs = self.out_dir / "concatenated_contigs.fasta"
        self.out_dir.mkdir(parents=True, exist_ok=True)

    #def tearDown(self):
    #    remove_folder(self.run_project_out_dir)

    def test_curation_same_species(self):
        args = {"assembly_fpath": self.contigs, "reference_input": self.reference_sequence,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        assembly_equivalencies = get_genomic_coordinates_from_assembly(nucmer_results["output_file"])
        assert assembly_equivalencies == [{'assembly_start': 1, 'assembly_end': 84170, 'ref_start': 1, 'ref_end': 84170, 'contig_id': 'NC_000932.1:1-84170', 'strand': '+'}, 
                                          {'assembly_start': 1, 'assembly_end': 26264, 'ref_start': 84171, 'ref_end': 110434, 'contig_id': 'NC_000932.1:84171-110434', 'strand': '+'}, 
                                          {'assembly_start': 1, 'assembly_end': 17780, 'ref_start': 110435, 'ref_end': 128214, 'contig_id': 'NC_000932.1:110435-128214', 'strand': '+'}, 
                                          {'assembly_start': 1, 'assembly_end': 26264, 'ref_start': 128215, 'ref_end': 154478, 'contig_id': 'NC_000932.1:84171-110434', 'strand': '-'}]
    
    def test_calculate_repeats_regions(self):
        args = {"assembly_fpath": self.reference_sequence, "reference_input": self.reference_sequence,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        inverted_repeats_regions = (calculate_repeats_regions(nucmer_results["output_file"]))

        args = {"assembly_fpath": self.overlapping_contigs, "reference_input": self.reference_sequence,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        biggest_inverted_repeat_contig = find_biggest_inverted_repeat_sequence(inverted_repeats_regions, nucmer_results["output_file"])
        print(biggest_inverted_repeat_contig)
        args = {"assembly_fpath": self.overlapping_contigs, "reference_input": self.overlapping_contigs,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        not_overalpping_coordinates = retrieve_not_overlapping_coordinates_from_contigs(biggest_inverted_repeat_contig, nucmer_results["output_file"])
        write_not_overlapping_regions(str(self.overlapping_contigs), self.not_overlapping_regions_seq, not_overalpping_coordinates)

    def test_curation_same_species_with_overlpas(self):
        args = {"assembly_fpath": self.overlapping_contigs, "reference_input": self.overlapping_contigs,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        assembly_equivalencies = get_genomic_coordinates_from_assembly(nucmer_results["output_file"])
        print(assembly_equivalencies)
        args = {"assembly_fpath": self.not_overlapping_regions_seq, "reference_input": self.reference_sequence,
                "out_dir": self.run_project_out_dir}
        nucmer_results = run_nucmer(args)
        assembly_equivalencies = get_genomic_coordinates_from_assembly(nucmer_results["output_file"])
        print(assembly_equivalencies)
        orient_and_concatenate_contigs(str(self.not_overlapping_regions_seq), str(self.concatenated_contigs), str(self.oriented_contigs), nucmer_results["output_file"])
        #assert assembly_equivalencies == [{'assembly_start': 885, 'assembly_end': 86794, 'ref_start': 1, 'ref_end': 85910, 'contig_id': 'NC_000932.1:1-84170', 'strand': '+'}, 
        #                                  {'assembly_start': 19521, 'assembly_end': 19520, 'ref_start': 84171, 'ref_end': 85910, 'contig_id': 'NC_000932.1:110435-128214', 'strand': '+'}, 
                                        #   {'assembly_start': 1741, 'assembly_end': 26264, 'ref_start': 84171, 'ref_end': 110434, 'contig_id': 'NC_000932.1:84171-110434', 'strand': '+'}, 
                                        #   {'assembly_start': 885, 'assembly_end': 885, 'ref_start': 109551, 'ref_end': 110435, 'contig_id': 'NC_000932.1:1-84170', 'strand': '+'}, 
                                        #   {'assembly_start': 1, 'assembly_end': 17784, 'ref_start': 110435, 'ref_end': 128218, 'contig_id': 'NC_000932.1:110435-128214', 'strand': '+'}, 
                                        #   {'assembly_start': 1, 'assembly_end': 880, 'ref_start': 128215, 'ref_end': 129098, 'contig_id': 'NC_000932.1:1-84170', 'strand': '-'}, 
                                        #   {'assembly_start': 1, 'assembly_end': 25380, 'ref_start': 128215, 'ref_end': 154478, 'contig_id': 'NC_000932.1:84171-110434', 'strand': '-'}, 
                                        #   {'assembly_start': 85055, 'assembly_end': 85054, 'ref_start': 152739, 'ref_end': 154478, 'contig_id': 'NC_000932.1:1-84170', 'strand': '-'}, 
                                        #   {'assembly_start': 17781, 'assembly_end': 17780, 'ref_start': 152739, 'ref_end': 154478, 'contig_id': 'NC_000932.1:110435-128214', 'strand': '-'}]
        #create_curated_genome(args, assembly_equivalencies)
            

        