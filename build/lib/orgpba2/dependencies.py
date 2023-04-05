from os import environ as env
from os.path import (exists, join)
from shutil import which


def check_binary_in_user_envs(user_env):
    #This function checks if user has defined a path
    #for program needed to run orgpba2
    program , _ = user_env.items()
    binaries = user_env["binaries"]
    user_path = user_env["user_path"]
    #Try to find if variable is user_path is defined
    try:
        user_path_check = env[user_path]
    except KeyError:
        msg = "{} is not defined"
        msg +="checking if {} is in $PATH"
        print(msg.format(program))
        return False
    #If user_path is defined, try to find all binaries
    #needed to run the program
    else:
        for binary in binaries:
           if not exists(join(user_path_check, binary)):
               msg = "{} doesn't contain"
               msg += " {} binary, checking if "
               msg += "{} is defined in $PATH"
               print(msg.format(user_path_check, binary,
                                program))
               return False
    return True

# def check_binary_requirements(binaries=BINARIES_REQUIREMENTS):
#     for dependency, features in binaries.items():

        

