import unittest

from pathlib import Path

from src.coverage import sort_sequences_by_position, get_reads_by_coverage


class TestCoverture(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.alignments_path = test_path / "01_mapping_minimap2.paf"

    def test_sort_sequences_by_position(self):
        sorted_alignments = sort_sequences_by_position(self.alignments_path)
        # print(sorted_alignments[0])
        # print(sorted_alignments[1])
        # print(sorted_alignments[-2])
        # print(sorted_alignments[-1])

    def test_retrieve_reads_by_coverture(self):
        reference_length = 154474
        sorted_alignments = sort_sequences_by_position(self.alignments_path)
        desired_coverage = 20
        fraction_window = 0.9
        reads_retrieved, coverages = get_reads_by_coverage(sorted_alignments, reference_length, desired_coverage, fraction_window=fraction_window)
        
    

