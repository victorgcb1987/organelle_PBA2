import shutil
import unittest

from pathlib import Path
from src.canu import run_canu
from src.seqs import get_seq_length, count_seqs, get_seqs_length
from src.utils import parse_executable_options

class TestCanu(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "data"
        self.test_seqs = self.test_path / "seqs_to_assemble.fq.gz"
        self.reference = self.test_path / "artha_chloro_ref.fasta.gz"

    def tearDown(self):
        shutil.rmtree(self.test_path / "03_canu")
    
    def test_canu(self):
        organule_size = get_seq_length(self.reference)
        canu_specific_arguments = parse_executable_options("-p=test")
        canu_arguments = {"genome_size": organule_size,
                          "canu_options": canu_specific_arguments,
                          "number_threads": 4,
                          "out_dir": self.test_path,
                          "sequence_technology": "pacbio",
                          "seqs_input": self.test_seqs}
        canu_results = run_canu(canu_arguments)
        output = canu_results["output_files"][0] / "test.contigs.fasta"
        assert count_seqs(output) == 2
        assert get_seqs_length(output) == 170021

if __name__ == "__main__":
    unittest.main()