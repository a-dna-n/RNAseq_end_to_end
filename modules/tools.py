import sys
import os
from collections import Counter
import subprocess
from typing import TypeAlias, Union #, Literal
import  shlex

File_or_Dir: TypeAlias = Union[str | bytes | os.PathLike]

# cat grep hostname wget which zcat
# abra2 minimap2 pysradb regtools samtools spades STAR 

def log_message(*args, exit_now: bool = False, **kwargs):
    print(*args, **kwargs, file=sys.stderr)
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
    splitlines: bool = True,
):
    
    with subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
        output = process.communicate()[0]
    output = output.decode("utf-8") #.rstrip(')
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except:
        log_message(f"Execution failed: {command}", exit_now = exit_on_error)
        return []
    output = output.decode("utf-8") #.rstrip(')
    """
    if splitlines:
        return output.splitlines()
    return output


def verify_columns_present_in_dataframe(
    *, data, columns: str | list[str], source: str
):
    # data = DataFrame
    # columns = string, set or list of strings
    # if any columns are missing, outputs error message and exits
    if isinstance(columns, str):
        columns = [columns]
    if missing := set(columns) - set(data.columns):
        log_message(f"Column(s) missing in {source}:", *missing, sep="\n\t", exit_now=True)

def verify_columns_absent_in_dataframe(
    *, data, columns: str | list[str], source: str
):
    # data = DataFrame
    # columns = string, set or list of strings
    # if any columns are missing, outputs error message and exits
    if isinstance(columns, str):
        columns = [columns]
    if conflict := set(columns) & set(data.columns):
        log_message(f"Column(s) found {source}:", *conflict, sep="\n\t", exit_now=True)

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
    if missing := [exe for exe in exes if not execute_command(f"which {exe}", exit_on_error=False)]:
        log_message("Missing exes:", *missing, sep="\n\t",  exit_now = exit_on_error)

    #for exe in exes:
    #    execute_command(f"which {exe}", exit_on_error=exit_on_error)

def test_libraries(libraries : str | list[str], *, exit_on_error: bool = True):
    """
    Test if libraries  are installed.

    Args:
        libraries : string or list of strings
        exit_on_error: stop if any package on the list is not installed.
    """
    from importlib.util import find_spec
    if missing := [p for p in libraries if not find_spec(p)]:
        log_message("Missing libraries:", *missing, sep="\n\t",  exit_now = exit_on_error)

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


def check_dir_write_access(paths: File_or_Dir | list[File_or_Dir]):
    # if dir exists, checks for write access
    # otherwise, tries to create it
    # returns dir or fails
    if isinstance(paths, File_or_Dir):
        paths = [paths]
    for p in paths:
        p = p.rstrip("/")
        if os.path.exists(p):
            if os.access(p, os.W_OK):
                log_message(f"{p} exists, write access OK")
            else:
                log_message(f"No write access to {p}", exit_now=True)
        else:
            try:
                os.makedirs(p)
                log_message(f"{p} created")
            except:
                log_message(f"Cannot create {p}", exit_now=True)
            

def check_if_dir_is_mounted(directory: File_or_Dir):
    for line in execute_command("mount"):
        temp = line.split(" ")
        if len(temp) > 1 and temp[1] == "on" and directory.startswith(temp[2]):
            log_message(f"{temp[2]} is mounted")
            return
    log_message(f"{directory} is not on a mounted drive", exit_now=True)


def is_a_gene_bam(bamfile: File_or_Dir):
    if ".realigned." in bamfile:
        log_message(f"Unexpected input: {bamfile} is a realigned BAM", exit_now=True)
    if ".genes." in bamfile:
        return True
    return False


def gene_bam(bamfile: File_or_Dir):
    # input = name of standard bam
    # output = corresponding gene-specific bam
    if is_a_gene_bam(bamfile):
        log_message(
            f"Unexpected input: {bamfile} is a region-specific BAM file", exit_now=True
        )
    return re.sub(r"\.bam$", ".genes.bam", bamfile)


def realigned_gene_bam(bamfile: File_or_Dir):
    # input = gene bam file
    # output = corresponding gene-specific bam
    if is_a_gene_bam(bamfile):
        return re.sub(".genes.bam", ".genes.realigned.bam", bamfile)
    log_message(
        f"Unexpected input: {bamfile} is not a region-specific BAM file.", exit_now=True
    )


def get_read_type_from_a_bamfile(bamfile: File_or_Dir): #, check_if_bamfile_exists: bool = True):
    # bamfile = string
    # checks flag of first alignment, returns "paired" or constants.read_type.single
    verify_that_paths_exist(bamfile)
    samtools_flag = execute_command(
        f"samtools view {bamfile} | head -n 1 | cut -f 2", exit_on_error=True
    )[0]
    if isinstance(samtools_flag, str):
        samtools_flag_translation = execute_command(
            f"samtools flags {samtools_flag}", exit_on_error=True
        )[0]
        if (
            "PAIRED" in samtools_flag_translation
        ):  # execute_command(f"samtools flags {samtools_flag} | grep PAIRED", exit_on_error=True):
            return constants.read_type.paired
        else:
            return constants.read_type.single
    log_message(
        f"Invalid output from samtools view {bamfile}:\n{samtools_flag}", exit_now=True
    )

def get_species_from_a_bam_file(bamfile: File_or_Dir):
    for line in execute_command(f"samtools view -H {bamfile}", exit_on_error=True):
        if m := re.search(r'.*?genomeDir\s+(.*?)\s', line):
            for s in constants.known_species:
                if s in m.groups()[0]:
                    return s
            log_message(f"{bamfile}: genomeDir but no valid species:\n" + line.rstrip())
            return None
    log_message(f"No species was found for {bamfile}")
    return None

def get_species_for_bamfiles(
    bamfiles: list[File_or_Dir],
):  # , known_species: constants.known_species):
    # input = bam files as list of strings
    # expects to find species in STAR genomeDir parameter from samtools header
    # @PG	ID:STAR	PN:STAR	VN:2.7.10b	CL:STAR   ... --genomeDir ${HOME}/star/human.GRCh38_109.75 ...
    # fails if no species is found or multiple species are found
    # otherwise returns a dict

    verify_that_paths_exist(bamfiles)
    bam_to_species = {bamfile: get_species_from_a_bam_file(bamfile) for bamfile in bamfiles}
    bam_to_species = {k:v for k in bam_to_species if k is not None}
    return bam_to_species


def SRA_simplify_exp_desc(x):
    # GSM4503604: NSC-CB660-TERT sgTP53 + sgCDKN2A + sgPTEN +sgNF1 CLONE 1 replicate 1; Homo sapiens; RNA-Seq
    desc = str(x["experiment_title"])
    for col in set(["organism_name", "library_strategy"]) & (set(x.keys())):
        temp = "; " + str(x[col])
        desc = re.sub(temp, "", desc)
    for col in set(["library_name", "experiment_alias"]) & (set(x.keys())):
        temp = str(x[col]) + ": "
        desc = re.sub(temp, "", desc)
    desc = re.sub(" ", "_", desc)
    return desc

def get_single_species_for_bamfiles(*, bamfiles: list[File_or_Dir]):  # , check_if_bamfiles_exist=True):
    # bamfiles = list of bam files
    # reports error and stops if there's more than one species for a list of BAM files
    # otherwise returns species as a string
    assert isinstance(bamfiles, list), "bamfiles is not a list"
    # if check_if_bamfiles_exist:
    verify_that_paths_exist(bamfiles)
    # species_for_bams = list(set([execute_command(f"{get_species_for_bam_sh} {bam}", exit_on_error=True) for bam in bamfiles]))
    species_for_bams = list(set(get_species_for_bamfiles(bamfiles=bamfiles).values()))
    if len(species_for_bams) == 1:
        return species_for_bams[0]
    else:
        log_message("Invalid number of species for bam files:", *species_for_bams, sep="\n", exit_now=True)


def get_read_type_from_bamfiles(*, bamfiles: list[File_or_Dir],
):  # , check_if_bamfiles_exist=True):
    # bamfiles = list of bam files
    # reports error and stops if there's more than one read type
    # otherwise returns read type (paired or single) as a string
    #assert isinstance(bamfiles, list), "bamfiles is not a list"
    # if check_if_bamfiles_exist:
    verify_that_paths_exist(bamfiles)
    # read_types = list(set([execute_command(f"{determine_read_type_sh} {bam}", exit_on_error=True) for bam in bamfiles]))
    read_types = list(
        set(
            [
                get_read_type_from_a_bamfile(bamfile) #=bam, check_if_bamfile_exists=False)
                for bam in bamfiles
            ]
        )
    )
    """
    if len(read_types) == 1:
        return read_types[0]
    else:
        temp = " ".join(read_types)
        log_message(f"Multiple read types : {temp}", exit_now=True)
    """

def make_gene_bam(*, inputbam: File_or_Dir, outputbam: File_or_Dir, regions: File_or_Dir, execute: bool = False):
    verify_that_paths_exist([inputbam, regions])
    exit_if_files_exist(outputbam)
    command = f"samtools view -o {outputbam} --region {regions} {inputbam} && samtools index {outputbam}"
    if execute:
        execute_command(command=command, log_command=True, exit_on_error=True)
        for x in [outputbam, f"{outputbam}.bai"]:
            if not os.path.exists(x):
                log_message(f"Status OK but {x} not found", exit_now=True)
    else:
        return command

def check_for_file_xor_stdin(inputfile: File_or_Dir | None):
    if inputfile is None:
        if detect_stdin():
            return
        log_message("Specify an input file or stdin.", exit_now=True)
    elif detect_stdin():
        log_message("Specify an input file or stdin, but not both.", exit_now=True)
        verify_that_paths_exist(inputfile)


def unique_row_headers(data: list[list]):
    # input list of lists
    # makes first column IDs unique by adding 1, 2 etc. if repeated
    # assert isinstance(data, list), "data is not a list"
    # assert all(isinstance(x, list) for x in data), "data is not a list of lists"

    ctr = Counter([row[0] for row in data])
    repeat = {x: 1 if y > 1 else 0 for (x, y) in ctr.items()}
    for i, row in enumerate(data):
        id = row[0]
        if repeat[id]:
            data[i][0] = f"{id}_{repeat[id]}"
            repeat[id] += 1
    return data
    """
    input:
        characteristics
        characteristics
        characteristics

    output:
        characteristics_1
        characteristics_2
        characteristics_3

    """

