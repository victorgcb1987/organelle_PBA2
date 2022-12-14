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
                      complement)
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
        self.sequence_to_check_origin = self.test_path / "artha_chloro_ref.fasta"
        self.circular_sequence = self.test_path / "circular_sequence.fasta"
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
       

    def test_if_fasta_or_fastq(self):
        #Check if file is a fasta or fastq
        assert sequence_kind(self.fastq_path) == "fastq"

        seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
        assert sequence_kind(self.fasta_path) == "fasta"

        seq_path = self.test_path / "dummy"
        self.assertRaises(RuntimeError, sequence_kind, seq_path)


    def test_seq_count(self):
        #checks number of seqs in file
        seq_path = self.test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
        assert count_seqs(seq_path) == 24903
        seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
        assert count_seqs(seq_path) == 1



    def test_get_seqs_ids_from_paf(self):
        #Check if identifiers from paf file are properly collected
        assert len(get_seqs_id_from_paf_file(self.paf_path)) == 1061


    def test_write_seqs_from_ids(self):
        #Check if reads from paf file are properly written
        seq_ids = get_seqs_id_from_paf_file(self.paf_path)
        write_seqs_from_seqs_id(seq_ids, self.fastq_path, self.seqs_out)
        assert count_seqs(self.seqs_out) == 1061
        self.seqs_out.unlink()


    def test_find_origin(self):
        #Check if circularity is found between assemblies
        blast_arguments = {"kind": "blastn",
                           "task": "megablast",
                           "outmft": "-outfmt '6 std qlen slen'" }
        options = {"out_dir": self.test_path}
        blastn_input = run_blast(self.blastdb_path, self.blastdb_path,
                                 options, blast_arguments=blast_arguments)
        assert blastn_input["return_code"] == 0

        origin = find_origin(blastn_input["output_file"], margin=10)
        assert origin["reference"] == 0
        assert origin["strand"] == "+"


        # Comparing origin between original assembly and reverse complement assembly
        reverse_comp_fhand = open(self.reverse_comp, "w")
        record = SeqIO.read(self.blastdb_path, "fasta")
        reverse_comp_fhand.write(record.id)
        reverse_comp_fhand.write(str(record.seq.reverse_complement()))
        reverse_comp_fhand.close()
        blastn_input = run_blast(self.reverse_comp, self.blastdb_path, 
                                 options, blast_arguments=blast_arguments)
        assert blastn_input["return_code"] == 0
        origin = find_origin(blastn_input["output_file"], margin=10)
        assert origin["reference"] == 0
        assert origin["query"] == 154480
        assert origin["strand"] == "-"


    def test_reverse_complement(self):
        #Check reverse complement
        seq = "ATCG"
        assert reverse_complement(seq) == "CGAT"


    def test_complement(self):
        #Check complement
        seq = "ATCG"
        assert complement(seq) == "TAGC"


    def test_circularity(self):
        assembly = self.circular_sequence
        options = {"out_dir": self.test_path}
        circularity = find_circularity(assembly, options, overlap_length=20)
        print(circularity)
        assert circularity["overlap_start"] == (1, 64)
        assert circularity["overlap_end"] == (125, 188)
        assert circularity["overlap_identity"] == 100
        self.assertFalse(find_circularity(assembly, options, overlap_length=100))

    

    # def test_remove_circularity_redundance(self):
    #     blast_arguments = {"kind": "blastn",
    #                        "task": "megablast",
    #                        "outmft": "-outfmt '6 std qlen slen'" }
    #     options = {"out_dir": self.test_path}
    #     blastn_input = run_blast(self.circular_sequence, self.sequence_to_check_origin,
    #                              options, blast_arguments=blast_arguments)
        
    #     origin = find_origin(blastn_input["output_file"], margin=10)
    #     circularity = find_circularity(self.circular_sequence, options, overlap_length=100)
    #     no_redundant_seq = remove_circularity_redundancy(self.circular_sequence, circularity)
    #     print(origin)
    #     print(circularity)
    #     reconstructed_sequence = reconstruct_assembly_from_origin(no_redundant_seq, origin)
    #     print(len(reconstructed_sequence))
    #     print(">reconstructed_seq")
    #     print(reconstructed_sequence)
    # #     assert reconstructed_sequence == "ATCGCATACGATCAGATCGCATATATATTATCGCTAGCTGACTATCGCAGCATCAGTCAccgtaagatagacgcacgtcccgtaagatagacgcacgtcccgtaagatagacgcacgtcCCCCG"


if __name__ == "__main__":
    unittest.main()