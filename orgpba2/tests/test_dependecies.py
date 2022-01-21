import unittest
from os import environ as env
from pathlib import Path

from orgpba2.config  import BINARIES_REQUIREMENTS
from orgpba2.dependencies import (check_binary_in_user_envs, 
                                  check_is_binary_is_in_PATH)

# This test unit checks if programs needed for runnning 
# orgpba are reachable in user defined enviromental 
# variables or in $PATH. User env path variables names
# checked for running pbaorg2 are in config.py

class TestDependencies(unittest.TestCase):

    def test_binary_in_user_envs(self):
        # Get dummy dummy user env path
        # test data folder has dummy files as binaries
        test_path = Path(__file__).parent.absolute()
        binary_path = test_path / "data"    
        env["DUMMY_PATH"] = str(binary_path)

        # Both dummy binary and user env are defined
        dummy_correct = {"dummy": 
                         {"binaries": ["dummy"],
                          "user_path": "DUMMY_PATH"}}
        assert check_binary_in_user_envs(dummy_correct), "NO"

        # user env path is defined but dummy binary doesn't exist in
        # in that path
        dummy_not_binary = {"dummy": 
                            {"binaries": ["not_found"],
                             "user_path": "DUMMY_PATH"}}
        self.assertFalse(check_binary_in_user_envs(dummy_not_binary))

        # binary is defined but user env path is not
        dummy_not_user_path = {"dummy": 
                               {"binaries": ["dummy"],
                                "user_path": "NOT_FOUND"}}
        self.assertFalse(check_binary_in_user_envs(dummy_not_user_path))

        # check when multiple binaries are needed for running a program
        dummy_multiple_binaries= {"dummy": 
                                  {"binaries": ["dummy", "dummy2"],
                                   "user_path": "DUMMY_PATH"}}
        assert check_binary_in_user_envs(dummy_multiple_binaries)

    
    def test_is_binary_is_in_PATH(self):
        binaries = ["cat", "grep"]
        assert check_is_binary_is_in_PATH(binaries)

        binary = ["not_found"]
        self.assertFalse(check_is_binary_is_in_PATH(binary))

if __name__ == "__main__":
    unittest.main()


