import unittest

from pathlib import Path

from orgpba2.seqs import (sequence_kind, count_seqs,
                          get_seqs_id_from_paf_file)

#sequence tools
class TestSeqUtils(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "data"

    def test_if_fasta_or_fastq(self):
        #Check if file is a fasta or fastq
        seq_path = self.test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
        assert sequence_kind(seq_path) == "fastq"

        seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
        assert sequence_kind(seq_path) == "fasta"

        seq_path = self.test_path / "dummy"
        self.assertRaises(RuntimeError, sequence_kind, seq_path)

    def test_seq_count(self):
        #checks number of seqs in file
        seq_path = self.test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
        assert count_seqs(seq_path) == 24903

        seq_path = self.test_path / "artha_chloro_ref.fasta.gz"
        assert count_seqs(seq_path) == 1

    def test_get_seqs_ids_from_paf(self):
        paf_path = self.test_path / "test.paf"
        assert len(get_seqs_id_from_paf_file(paf_path)) == 1061






if __name__ == "__main__":
    unittest.main()