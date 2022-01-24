import unittest

from pathlib import Path

from orgpba2.seqs import (sequence_kind, count_seqs)

#sequence tools
class TestSeqUtils(unittest.TestCase):

    def test_if_fasta_or_fastq(self):
        #Check if file is a fasta 
        test_path = Path(__file__).parent.absolute()

        seq_path = test_path / "data" / "artha_pacbioSRR1284093_c025k.fastq.gz"
        assert sequence_kind(seq_path) == "fastq"

        seq_path = test_path / "data" /  "artha_chloro_ref.fasta.gz"
        assert sequence_kind(seq_path) == "fasta"

        seq_path = test_path / "data" / "dummy"
        self.assertRaises(RuntimeError, sequence_kind, seq_path)

    def test_seq_count(self):
        #checks number of seqs in file
        test_path = Path(__file__).parent.absolute()

        seq_path = test_path / "data" / "artha_pacbioSRR1284093_c025k.fastq.gz"
        assert count_seqs(seq_path) == 24903

        seq_path = test_path / "data" /  "artha_chloro_ref.fasta.gz"
        assert count_seqs(seq_path) == 1





if __name__ == "__main__":
    unittest.main()