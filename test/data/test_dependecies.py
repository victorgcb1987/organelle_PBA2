import unittest

from os import environ as env
from os import join
from os.path import exists
from shutil import which

from organlle_PBA2.src.config import (BINARIES, ENV_VARIABLES)


class TestDependencies(unittest.TestCase):

    def test_minimap2_dependency(self):
        try:
            #Checking if user has defined $MINIMAP2_PATH
            minimap2_path = env[ENV_VARIABLES["minimap2"]]
            
        except KeyError:
            #If is not defined, wil try to find it in $PATH
            msg = "$MINIMAP2_PATH is not defined"
            msg +=" checking if minimap2 is in $PATH"
            print(msg)

        else:
            #If $MINIMAP2 is defined, 
            #try to find minimap2 binary on it
            msg = "$MINIMAP2_PATH2 doesn't contain"
            msg += " minimap2 binary, checking if "
            msg += "minimap2 is defined in $PATH"
            assert exists(join(minimap2_path, 
                               BINARIES["minimap2"])), msg

        #Else will try to find minimap2 binary in $PATH
        msg = "minimap2 is not defined in $PATH"
        msg += " minimap2 path must be either in $PATH"
        msg += " or $MINIMAP2_PATH"
        assert which(BINARIES["minimap2"]), msg

if __name__ == "__main__":
    unittest.main()


