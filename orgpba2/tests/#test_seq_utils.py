import unittest

from pathlib import Path
from shutil import rmtree as remove_folder

from Bio import SeqIO

from src.seqs import (sequence_kind, count_seqs,
                      get_seqs_id_from_paf_file,
                      write_seqs_from_seqs_id,
                      run_blast, find_origin,
                      find_circularity,
                      remove_circularity_redundancy,
                      reconstruct_assembly_from_origin,
                      reverse_complement,
                      find_blocks_breakpoints,
                      create_haplotype_breakpoints_sequences,
                      write_haplotypes_breakpoints)
from src.utils import folder_exists

#sequence tools
class TestSeqUtils(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "data"
        self.paf_path = self.test_path / "test.paf"
        self.fastq_path = self.test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
        self.fasta_path = self.test_path / "artha_chloro_ref.fasta.gz"
        self.blastdb_path = self.test_path / "artha_chloro_ref.fasta"
        self.seqs_out = self.test_path / "out_seqs.fasta"
        self.reverse_comp = self.test_path / "reverse_comp.fasta"
        self.assembly = self.test_path / "scaffolds.fasta"
        #self.sequence_to_check_origin = self.test_path / "original.fasta"
        #self.circular_sequence = self.test_path / "circular.fasta"
        self.sequence_to_check_origin = self.test_path / "artha_chloro_ref.fasta"
        self.circular_sequence = self.test_path / "scaffolds.fasta"
        self.haplotypes = self.test_path / "haplotypes.fasta"
        self.colinear_irs = self.test_path / "colinear_IRs.fa"

    # def tearDown(self):
    #     blast_dir = self.test_path / "00_blasts"
    #     if folder_exists(blast_dir):
    #         remove_folder(blast_dir)
    #     for _file in self.test_path.glob("*.fasta.*"):
    #         if not str(_file).endswith(".gz"):
    #             file_path = self.test_path / _file
    #             file_path.unlink()
    #     blast_dir = self.test_path / "00_blasts"
       


    # def test_if_fasta_or_fastq(self):
    #     #Check if file is a fasta or fastq
    #     assert sequence_kind(self.fastq_path) == "fastq"

    #     seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
    #     assert sequence_kind(self.fasta_path) == "fasta"

    #     seq_path = self.test_path / "dummy"
    #     self.assertRaises(RuntimeError, sequence_kind, seq_path)

    # def test_seq_count(self):
    #     #checks number of seqs in file
    #     seq_path = self.test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
    #     assert count_seqs(seq_path) == 24903
    #     seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
    #     assert count_seqs(seq_path) == 1

    # def test_get_seqs_ids_from_paf(self):
    #     assert len(get_seqs_id_from_paf_file(self.paf_path)) == 1061

    # def test_write_seqs_from_ids(self):
    #     seq_ids = get_seqs_id_from_paf_file(self.paf_path)
    #     write_seqs_from_seqs_id(seq_ids, self.fastq_path, self.seqs_out)
    #     assert count_seqs(self.seqs_out) == 1061
    #     self.seqs_out.unlink()

    # def test_find_origin(self):
    #     blast_arguments = {"kind": "blastn",
    #                        "task": "megablast",
    #                        "outmft": "-outfmt '6 std qlen slen'" }
    #     options = {"out_dir": self.test_path}
    #     blastn_input = run_blast(self.blastdb_path, self.blastdb_path,
    #                              options, blast_arguments=blast_arguments)
    #     assert blastn_input["return_code"] == 0

    #     origin = find_origin(blastn_input["output_file"], margin=10)
    #     assert origin["reference"] == 1
    #     assert origin["strand"] == "+"


    #     reverse_comp_fhand = open(self.reverse_comp, "w")
    #     record = SeqIO.read(self.blastdb_path, "fasta")
    #     reverse_comp_fhand.write(record.id)
    #     reverse_comp_fhand.write(str(record.seq.reverse_complement()))
    #     reverse_comp_fhand.close()
        
    #     blastn_input = run_blast(self.reverse_comp, self.blastdb_path, 
    #                              options, blast_arguments=blast_arguments)
    #     assert blastn_input["return_code"] == 0
    #     origin = find_origin(blastn_input["output_file"], margin=10)
    #     assert origin["reference"] == 1
    #     assert origin["strand"] == "-"

    # def test_circularity(self):
    #     assembly = self.circular_sequence
    #     options = {"out_dir": self.test_path}
    #     circularity = find_circularity(assembly, options, overlap_length=60)
    #     assert circularity["overlap_start"] == (1, 64)
    #     assert circularity["overlap_end"] == (125, 188)
    #     assert circularity["overlap_identity"] == 100
    #     self.assertFalse(find_circularity(assembly, options, overlap_length=100))

    # def test_reverse_complement(self):
    #     seq = "ATCG"
    #     assert reverse_complement(seq) == "CGAT"

    def test_remove_circularity_redundance(self):
        blast_arguments = {"kind": "blastn",
                           "task": "megablast",
                           "outmft": "-outfmt '6 std qlen slen'" }
        options = {"out_dir": self.test_path}
        blastn_input = run_blast(self.circular_sequence, self.sequence_to_check_origin,
                                 options, blast_arguments=blast_arguments)
        
        origin = find_origin(blastn_input["output_file"], margin=10)
        circularity = find_circularity(self.circular_sequence, options, overlap_length=100)
        no_redundant_seq = remove_circularity_redundancy(self.circular_sequence, circularity)
        print(origin)
        print(circularity)
        reconstructed_sequence = reconstruct_assembly_from_origin(no_redundant_seq, origin)
        print(len(reconstructed_sequence))
        print(">reconstructed_seq")
        print(reconstructed_sequence)
    #     assert reconstructed_sequence == "ATCGCATACGATCAGATCGCATATATATTATCGCTAGCTGACTATCGCAGCATCAGTCAccgtaagatagacgcacgtcccgtaagatagacgcacgtcccgtaagatagacgcacgtcCCCCG"

    def test_heteroplasmy(self):
        options = {"out_dir": self.test_path}
        breakpoints = find_blocks_breakpoints(self.blastdb_path, options)
        # assert breakpoints["LSC_IR"] == 84171
        # assert breakpoints["IR_SSC"] == 110434
        # assert breakpoints["SSC_IR"] == 128215
        # assert breakpoints["IR_LSC"] == 154478
        # breakpoints = find_blocks_breakpoints(self.colinear_irs, options)
        # assert breakpoints["LSC_IR"] == 84171
        # assert breakpoints["IR_SSC"] == 110434
        # assert breakpoints["SSC_IR"] == 128215
        # assert breakpoints["IR_LSC"] == 154478
        # out_fhand = open("Andropogon_gerardii.blocks.fasta", "w")
        # sequence = open(self.blastdb_path)
        # sequence.readline()
        # sequence = "".join([line.rstrip() for line in sequence])
        # # breakpoints = {"LSC_IRa": 4, "IRa_SSC": 8,  "SSC_IRb": 12, "IRb_LSC": 20}
        # # sequence =  "ATCATTGCATACGATGGACG"
        # parts = create_haplotype_breakpoints_sequences(sequence, breakpoints)
        # write_haplotypes_breakpoints(parts, breakpoints, out_fhand)




if __name__ == "__main__":
    unittest.main()