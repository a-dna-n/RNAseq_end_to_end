import sys
import os
import re
from collections import OrderedDict, Counter
import subprocess
from typing import TypeAlias, Union #, Literal

FileName: TypeAlias = Union[str | bytes | os.PathLike]

# cat grep hostname wget which zcat
# abra2 minimap2 pysradb regtools samtools spades STAR 

def log_message(message, *, exit_now: bool = False):
    print(message, file=sys.stderr)
    if exit_now:
        sys.exit(1)


def detect_stdin():
    if sys.stdin.isatty() == 0:
        return True
    else:
        return False

def _caller():
    import inspect
    return inspect.currentframe().f_code.co_name

def execute_command(
    command: str,
    *,
    exit_on_error: bool = False,
    log_command: bool = False,
    suppress_error_messages: bool = False,
):
    if log_command:
        log_message(f"# Executing {command}")
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except:
        if not suppress_error_messages:
            log_message(f"Execution failed: {command}")
        if exit_on_error:
            sys.exit(1)
        return []
    return [line.rstrip().decode("CP437") for line in output.splitlines()]


def verify_columns_present_in_dataframe(
    *, data, columns: str | list[str], source: str | None
):
    # data = DataFrame
    # columns = string, set or list of strings
    # if any columns are missing, outputs error message and exits
    if isinstance(columns, str):
        columns = [columns]
    if missing := columns - set(data.columns):
        missing = "\n".join(list(missing))
        if source:
            log_message(f"column(s) missing in {source}:")
        else:
            log_message("column(s) missing:")
        log_message(missing, exit_now=True)


def system_memory():
    return int(execute_command("free --mega | grep ^Mem")[0].split(" ")[-1])

def test_executables(exes: str | list[str], *, exit_on_error: bool = False):
    """
    Test if programs/scripts are executable.

    Args:
        exes: program (string) or programs (list of strings)
        exit_on_error: stop if any program on the list is not executable. Defaults to False.
    """
    # exes = string or list of strings
    # check if executable
    if isinstance(exes, str):
        exes = [exes]
    for exe in exes:
        execute_command(f"which {exe}", exit_on_error=exit_on_error)


def verify_that_paths_exist(path_or_paths: str | list[str], exit_on_error: bool = True):
    # checks if files and/or folders exist
    # returns true if all exist
    if isinstance(path_or_paths, str):
        path_list = [path_or_paths]
    else:
        path_list = path_or_paths
    if missing := [x for x in path_list if not os.path.exists(x)]:
        missing = "\n".join(missing)
        log_message(f"paths not found:\n{missing}", exit_now=exit_on_error)
        return False
    return True


def exit_if_files_exist(file_or_files: str | list[str]):
    # checks if files and/or folders exist
    # reports error and exits if any file exists
    if isinstance(file_or_files, str):
        file_list = [file_or_files]
    else:
        file_list = file_or_files
    if files_found := [x for x in file_list if os.path.exists(x)]:
        files_found = "\n".join(files_found)
        log_message(f"output file(s) found:\n{files_found}", exit_now=True)


def check_dir(dir=None):
    # if dir exists, checks for write access
    # otherwise, tries to create it
    # returns dir or fails
    assert isinstance(dir, str), "dir is not a string"
    dir = dir.rstrip("/")
    if os.path.exists(dir):
        if os.access(dir, os.W_OK):
            log_message(f"{dir} exists, write access OK")
            return dir
        log_message(f"No write access to {dir}")
    else:
        try:
            os.makedirs(dir)
            log_message(f"{dir} created")
            return dir
        except:
            log_message(f"Cannot create {dir}")
    for line in execute_command("mount"):
        temp = line.split(" ")
        if len(temp) > 1 and temp[1] == "on" and dir.startswith(temp[2]):
            log_message(f"{temp[2]} is mounted", exit_now=True)
    log_message(f"{dir} is not on a mounted drive", exit_now=True)


