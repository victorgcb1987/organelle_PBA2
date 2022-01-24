import unittest

from pathlib import Path

from orgpba2.utils import (file_exists, parse_args, file_is_compressed)

#Checks I/O tools
class TestUtils(unittest.TestCase):
    

    def test_parse_args(self):
        #parse arguments provided by the user at orgpba2 init
        #arguments must have the following syntax : arg1=value1,arg2=value3
        
        args = "a=1,b=2,c=3"
        assert parse_args(args) == {"a": "1", "b": "2", "c": "3"}
        
        args  = "a=1b=1;c=3"
        self.assertRaises(RuntimeError, parse_args, args)

    def test_file_exists(self):
        test_path = Path(__file__).parent.absolute()
        path = str(test_path / "data")
    
        assert file_exists(path, "dummy")

        self.assertFalse(file_exists(path, "notfound.exe"))

    def test_if_file_is_compressed(self):
        test_path = Path(__file__).parent.absolute()
        path = test_path / "data" / "artha_chloro_ref.fasta.gz"

        assert file_is_compressed(path)

        path = test_path / "data" / "dummy"

        self.assertFalse(file_is_compressed(path))







        

