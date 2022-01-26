from os import environ as env
from pathlib import Path
from shutil import which


def check_executable_in_user_envs(user_env, program="default"):
    #This function checks if user has defined a path
    #for program needed to run orgpba2
    executables = user_env["executables"]
    user_path = user_env["user_path"]

    #Try to find if variable is user_path is defined
    try:
        user_path_check = Path(env[user_path])
    except KeyError:
        msg = "user env {} is not defined for {}"
        print(msg.format(user_path, program))
        return False
    #If user_path is defined, try to find all executables
    #needed to run the program
    else:
        for executable in executables:
            full_path = Path(user_path_check) / executable
            check = full_path.absolute().exists()
            if not check:
               msg = "{} doesn't contain {} executable"
               print(msg.format(user_path_check, executable))
               return False
            else:
                print("{} found at {}".format(executable, full_path))
                return user_path_check

def check_is_executable_is_in_PATH(executables):
    msg = "{} is not defined in $PATH."
    msg += " {} path must be set either in $PATH"
    msg += " or set in user env path"
    for executable in executables:
        check = which(executable)
        if check:
            print("{} found at {}".format(executable, check))
        else:
            print(msg.format(executable, executable))
            return False
    return True

