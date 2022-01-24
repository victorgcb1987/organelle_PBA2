import unittest
from os import environ as env
from pathlib import Path

from orgpba2.config  import EXECUTABLES_REQUIREMENTS
from orgpba2.dependencies import (check_executable_in_user_envs, 
                                  check_is_executable_is_in_PATH)

# This test unit checks if executables needed for runnning 
# orgpba are reachable in user defined enviromental 
# variables or in $PATH. User env path variables names
# checked for running pbaorg2 are in config.py

class TestDependencies(unittest.TestCase):

    def test_executable_in_user_envs(self):
        # Get dummy dummy user env path
        # test data folder has dummy files as executables
        test_path = Path(__file__).parent.absolute()
        executable_path = test_path / "data"    
        env["DUMMY_PATH"] = str(executable_path)

        # Both dummy executable and user env are defined
        dummy_correct = {"dummy": 
                         {"executables": ["dummy"],
                          "user_path": "DUMMY_PATH"}}
        assert check_executable_in_user_envs(dummy_correct), "NO"

        # user env path is defined but dummy executable doesn't exist in
        # in that path
        dummy_not_executable = {"dummy": 
                            {"executables": ["not_found"],
                             "user_path": "DUMMY_PATH"}}
        self.assertFalse(check_executable_in_user_envs(dummy_not_executable))

        # executable is defined but user env path is not
        dummy_not_user_path = {"dummy": 
                               {"executables": ["dummy"],
                                "user_path": "NOT_FOUND"}}
        self.assertFalse(check_executable_in_user_envs(dummy_not_user_path))

        # check when multiple executables are needed for running a program
        dummy_multiple_executables= {"dummy": 
                                  {"executables": ["dummy", "dummy2"],
                                   "user_path": "DUMMY_PATH"}}
        assert check_executable_in_user_envs(dummy_multiple_executables)

    #Check if executable is in $PATH
    def test_is_executable_is_in_PATH(self):
        executables = ["cat", "grep"]
        assert check_is_executable_is_in_PATH(executables)

        executable = ["not_found"]
        self.assertFalse(check_is_executable_is_in_PATH(executable))    

if __name__ == "__main__":
    unittest.main()


