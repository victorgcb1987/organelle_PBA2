from os import environ as env
from os.path import (exists, join)
from shutil import which


def check_binary_in_user_envs(user_env):
    #This function checks if user has defined a path
    #for program needed to run orgpba2
    program = list(user_env.keys())[0]
    binaries = user_env[program]["binaries"]
    user_path = user_env[program]["user_path"]

    #Try to find if variable is user_path is defined
    try:
        user_path_check = env[user_path]
    except KeyError:
        msg = "user env {} is not defined for {}"
        print(msg.format(user_path, program))
        return False
    #If user_path is defined, try to find all binaries
    #needed to run the program
    else:
        for binary in binaries:
            full_path = join(user_path_check, binary)
            check = exists(full_path)
            if not check:
               msg = "{} doesn't contain {} binary"
               print(msg.format(user_path_check, binary))
               return False
            else:
                print("{} found at {}".format(binary, full_path))
    return True

def check_is_binary_is_in_PATH(binaries):
    msg = "{} is not defined in $PATH."
    msg += " {} path must be set either in $PATH"
    msg += " or set in user env path"
    for binary in binaries:
        check = which(binary)
        if check:
            print("{} found at {}".format(binary, check))
        else:
            print(msg.format(binary, binary))
            return False
    return True

# def check_binary_requirements(binaries=BINARIES_REQUIREMENTS):
#     for dependency, features in binaries.items():

        

