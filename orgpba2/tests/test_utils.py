import unittest

from pathlib import Path

from src.utils import (file_exists, 
                       parse_executable_options,
                       file_is_compressed,
                       _compress_file)

#Checks utils
class TestUtils(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "data"

    def test_file_exists(self):
        filepath = self.test_path / "dummy"
        assert file_exists(filepath)

        filepath = self.test_path / "notfound"
        self.assertFalse(file_exists(filepath))

    def test_if_file_is_compressed(self):
        filepath = self.test_path / "artha_chloro_ref.fasta.gz"
        assert file_is_compressed(filepath)

        filepath = self.test_path / "dummy"
        self.assertFalse(file_is_compressed(filepath))
    
    def test_file_compression(self):
        uncompressed_fpath = self.test_path / "uncompressed_test"
        compressed_fpath = self.test_path / "uncompressed_test.gz"
        _compress_file(uncompressed_fpath, compressed_fpath)
        assert compressed_fpath.stat().st_size < uncompressed_fpath.stat().st_size
        assert file_is_compressed(compressed_fpath)
        compressed_fpath.unlink()

class TestArgs(unittest.TestCase):

    def test_parse_executable_options(self):
        #parse arguments provided by the user at orgpba2 init
        #arguments must have the following syntax : arg1=value1,arg2=value3
        
        args = "a=1,b=2,c=3"
        assert parse_executable_options(args) == {"a": "1", "b": "2", "c": "3"}
        
        args  = "a=1b=1;c=3"
        self.assertRaises(RuntimeError, parse_executable_options, args)

if __name__ == "__main__":
    unittest.main()     
        
        

    
        


    







        

