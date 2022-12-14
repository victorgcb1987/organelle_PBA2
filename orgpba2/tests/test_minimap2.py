import unittest

from shutil import rmtree as remove_folder
from pathlib import Path

from src.minimap2 import run_minimap2

class TestMinimap2(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.reference_sequence = test_path / "artha_chloro_ref.fasta.gz"
        self.incorrect_reference = "incorrect.fasta"
        self.pacbio_sequences = test_path / "artha_pacbioSRR1284093_c025k.fastq.gz"
        self.run_project_out_dir = test_path / "test_out"
        self.out_minimap2_out = self.run_project_out_dir / "01_mapping_minimap2"
        self.out_minimap2_out.mkdir(parents=True, exist_ok=True)


    def tearDown(self):
        remove_folder(self.run_project_out_dir)


    def test_minimap2(self):
        # Checks that error messages are properly captured
        arguments = {"reference_input": self.incorrect_reference,
                     "sequences": self.pacbio_sequences,
                     "minimap2_options": {"-z": "incorrectvalue"},
                     "out_dir": self.run_project_out_dir,
                     "number_threads": 4,
                     "sequence_technology": "pacbio"}
        mininimap2_run = run_minimap2(arguments)
        assert mininimap2_run["return_code"] == 1
        assert "ERROR" in mininimap2_run["log_messages"] 

        #Check returncode and output fpath when minimap2 runs correctly
        arguments = {"reference_input": self.reference_sequence,
                     "sequences": self.pacbio_sequences,
                     "out_dir": self.run_project_out_dir,
                     "number_threads": 4,
                     "minimap2_options": {},
                     "sequence_technology": "pacbio"}
        mininimap2_run = run_minimap2(arguments)
        assert mininimap2_run["return_code"] == 0
        assert mininimap2_run["output_file"] == self.out_minimap2_out / "01_mapping_minimap2.paf"

        #Checks that module returns "already exists" if minimap2 output exists already
        mininimap2_run = run_minimap2(arguments)
        assert mininimap2_run["return_code"] == 0
        assert mininimap2_run["output_file"] == self.out_minimap2_out / "01_mapping_minimap2.paf"
        assert "already exists" in mininimap2_run["log_messages"]

        #Checks proper mapping output file extension if format i changed to sam
        arguments = {"reference_input": self.reference_sequence,
                     "sequences": self.pacbio_sequences,
                     "minimap2_options": {"-a": " "},
                     "out_dir": self.run_project_out_dir,
                     "number_threads": 4,
                     "sequence_technology": "pacbio"}
        mininimap2_run = run_minimap2(arguments)
        assert mininimap2_run["return_code"] == 0
        assert mininimap2_run["output_file"] == self.out_minimap2_out / "01_mapping_minimap2.sam"

        mininimap2_run = run_minimap2(arguments)
        assert mininimap2_run["return_code"] == 0
        assert mininimap2_run["output_file"] == self.out_minimap2_out / "01_mapping_minimap2.sam"
        assert "already exists" in mininimap2_run["log_messages"]
        

if __name__ == "__main__":
    unittest.main()