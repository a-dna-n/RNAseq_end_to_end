#!/usr/bin/env python

# /home/adnan/star/ref/ens.to.gene.mouse.if.unique

# functions for RNA-seq data - retrieval, alignment etc.

import sys
import os
import argparse
import re
import glob
from collections import OrderedDict, Counter
from textwrap import dedent
import subprocess
from types import SimpleNamespace
from typing import TypeAlias, Union #, Literal
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
#from modules.tools import *
from modules.arg_def import *

def is_a_gene_bam(bamfile: FileName):
    if ".realigned." in bamfile:
        log_message(f"Unexpected input: {bamfile} is a realigned BAM", exit_now=True)
    if ".genes." in bamfile:
        return True
    return False


def gene_bam(bamfile: FileName):
    # input = name of standard bam
    # output = corresponding gene-specific bam
    if is_a_gene_bam(bamfile):
        log_message(
            f"Unexpected input: {bamfile} is a region-specific BAM file", exit_now=True
        )
    return re.sub(r"\.bam$", ".genes.bam", bamfile)


def realigned_gene_bam(bamfile: FileName):
    # input = gene bam file
    # output = corresponding gene-specific bam
    if is_a_gene_bam(bamfile):
        return re.sub(".genes.bam", ".genes.realigned.bam", bamfile)
    log_message(
        f"Unexpected input: {bamfile} is not a region-specific BAM file.", exit_now=True
    )


def get_read_type_for_bamfile(bamfile: FileName, check_if_bamfile_exists: bool = True):
    # bamfile = string
    # checks flag of first alignment, returns "paired" or "single"
    if check_if_bamfile_exists:
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
            return "paired"
        else:
            return "single"
    log_message(
        f"Invalid output from samtools view {bamfile}:\n{samtools_flag}", exit_now=True
    )

def get_species_for_a_bam_file(bamfile: FileName):
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
    bamfiles: list[FileName],
):  # , known_species: constants.known_species):
    # input = bam files as list of strings
    # expects to find species in STAR genomeDir parameter from samtools header
    # @PG	ID:STAR	PN:STAR	VN:2.7.10b	CL:STAR   ... --genomeDir ${HOME}/star/human.GRCh38_109.75 ...
    # fails if no species is found or multiple species are found
    # otherwise returns a dict

    verify_that_paths_exist(bamfiles)
    species_for_bam = {}
    for bamfile in bamfiles:
        if species_found := get_species_for_a_bam_file(bamfile):
            species_for_bam[bamfile] = species_found
    return species_for_bam

r"""
species_found = False
#query = f"samtools view -H {bamfile} | grep genomeDir -m 1"
for line in execute_command("samtools view -H {bamfile}", exit_on_error=True):
    if m := re.search(r'.*?genomeDir\s+(.*?)\s', line):
        if m.groups()[0] in constants.known_species:
            species_for_bam[bamfile] = m.groups()[0]
            species_found = True
            break
if isinstance(temp, str):
    genomeDir = re.sub(r"", "", temp).split(" ")[0]
    found = [s for s in constants.known_species if s in genomeDir]
    if len(found) == 0:
        log_message(
            f"No recognized species in genomeDir {genomeDir} for BAM file {bamfile} in header line {temp}",
            exit_now=True,
        )
    elif len(found) > 1:
        found = ", ".join(found)
        log_message(
            f"Multiple species recognized in genomeDir {genomeDir} in BAM header line {temp}: {found}",
            exit_now=True,
        )
    else:
        species_for_bam[bamfile] = found[0]
        #for assembly, species in species_by_assembly.items():
        #    if assembly in genomeDir:
        #        return species
else:
    log_message(f"Invalid output from {query}", exit_now=True)

def sra_parse_exp(x):
    temp = str(x).split(": ")
    if len(temp) == 1: return temp[0]
    return temp[1].split(";")[0]
"""

def get_sra_metadata(*, id: str, delete_temp_output: bool = False):
    # get metadata for a single ID from sra via pysradb
    # returns dataframe
    
    db = SRAweb()
    return db.sra_metadata(id, detailed=True)
    """
    tempoutputfile = f"pysradb.{id}.temp"
    query = f"pysradb metadata --detailed {id} > {tempoutputfile}"
    df = db.sra_metadata("SRP016501")
    if os.path.exists(tempoutputfile):
        log_message(f"Using {tempoutputfile}")
    else:
        temp = execute_command(query)
        if not (os.path.exists(tempoutputfile)):
            log_message(f"No output from {query}")
            return None
    try:
        data = pd.read_csv(tempoutputfile, sep="\t", header=0, dtype=object)
    except:
        log_message(f"Error reading {tempoutputfile}. Query was: {query}")
        return None
    if delete_temp_output:
        os.remove(tempoutputfile)
        log_message(f"Deleted {tempoutputfile}.")
    return data
    """

def parse_ena_fastq_metadata(*, data):
    if not isinstance(data, pd.DataFrame):
        log_message(f"data is a {type(data)} instead of a dataframe")
        return None
    columns_expected = ["run_accession", "fastq_ftp"]
    if list(data.columns) != columns_expected:
        obs = ", ".join(data.columns)
        exp = ", ".join(columns_expected)
        log_message(f"Invalid columns: {obs}\nExpected: {exp}")
        return None
    data.columns = ["run_accession", "ena_fastq_ftp"]
    mask = data["ena_fastq_ftp"].isna()
    data = data[~mask].copy()
    data["ena_fastq_ftp_1"] = ""
    data["ena_fastq_ftp_2"] = ""
    data = data.set_index("run_accession")
    for srr in data.index:
        temp = data.at[srr, "ena_fastq_ftp"].split(";")
        if len(temp) > 2:
            r = data.at[srr, "ena_fastq_ftp"]
            log_message(f"\nWarning: unknown number of files in {r}\n")
        if len(temp) == 2:
            data.at[srr, "ena_fastq_ftp_1"] = temp[0]
            data.at[srr, "ena_fastq_ftp_2"] = temp[1]
            data.at[srr, "ena_fastq_ftp"] = ""
    drop = []
    for col in data.columns:
        if all(x == "" for x in data[col]):
            drop.append(col)
        else:
            data[col] = data[col].str.replace("ftp.sra.ebi.ac.uk/", "", regex=False)
    if drop:
        data.drop(columns=drop, inplace=True)
    return data.reset_index()


def get_enafastq_list(study_ID: str, *, delete_temp_output: bool = False):
    # gets list of ENA fastqs for ID
    # returns dataframe
    tempoutputfile = f"ena.{study_ID}.temp"
    query = f'wget -O {tempoutputfile} "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={study_ID}&result=read_run&fields=run_accession,fastq_ftp"'
    # https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP247494&result=read_run&fields=run_accession,fastq_ftp
    if os.path.exists(tempoutputfile):
        log_message(f"Using {tempoutputfile}")
    else:
        execute_command(query)
        if not (os.path.exists(tempoutputfile)):
            log_message(f"No output from {query} for {study_ID}")
            return None
    try:
        data = pd.read_csv(tempoutputfile, sep="\t", header=0, dtype=object)
    except:
        log_message(f"Error reading {tempoutputfile}. Query was: {query}")
        return None
    data = parse_ena_fastq_metadata(data=data)
    if data is None:
        log_message(f"Error parsing {tempoutputfile}. Query was: {query}")
        return None
    if delete_temp_output:
        os.remove(tempoutputfile)
        log_message(f"Deleted {tempoutputfile}.")
    return data

    """
    run_accession	fastq_ftp
    SRR12879058	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/058/SRR12879058/SRR12879058_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/058/SRR12879058/SRR12879058_2.fastq.gz
    SRR12879060	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/060/SRR12879060/SRR12879060.fastq.gz
    SRR12879061	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/061/SRR12879061/SRR12879061_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/061/SRR12879061/SRR12879061_2.fastq.gz
    SRR12879063	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/063/SRR12879063/SRR12879063.fastq.gz
    SRR12879064	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/064/SRR12879064/SRR12879064.fastq.gz
    SRR12879065	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/065/SRR12879065/SRR12879065.fastq.gz
    SRR12879066	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/066/SRR12879066/SRR12879066.fastq.gz
    SRR12879067	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/067/SRR12879067/SRR12879067.fastq.gz
    """


# def pysradb_conversions():
#    test_executables(exes=["pysradb"], exit_on_error=True)

def pysradb_parse_conversion(data, ID):
    if data is None:
        log_message(f"empty pysradb result for {ID}")
        return ""
    if not isinstance(data, pd.DataFrame):
        log_message(f"unexpected type {type(data)} for {ID}:\n{repr(data)}")
        return ""
    return data[data.columns[-1]].tolist()[0]

def pysradb_convert_ID(*, fn: str, id: str):
    """
    Calls pysradb to convert IDs using a specific
    fix order:
    > pysradb gsm-to-gse GSM4679562
    study_alias	study_accession
    GSE154783	SRP272683
    """
    cmd = f"pysradb {fn} {id}"
    log_message(cmd)
    temp = execute_command(cmd)

    prefix = fn[-3:].upper()
    # print(f"fn {fn} id {id} prefix {prefix}")
    if len(temp) == 2:
        new_ids = temp[-1].rstrip().split()
        # new_id = temp[-1].rstrip().split()[1]
        for new_id in new_ids:
            if new_id.startswith(prefix):
                log_message(f"{id} => {new_id}")
                return new_id
    temp = "\n".join(temp)
    log_message(f"Invalid output:\n{temp}")
    return None


def transpose_nested_list(*, data: list):
    # M x N list becomes N x M
    return [list(x) for x in zip(*data)]


def get_matrix_files_for_geo_series(geo_id: str, *, retrieve_matrix_files: bool = True):
    if not geo_id:
        log_message("Species a series  ID e.g. GSE154891", exit_now=True)
    # input = GEO series ID, e.g. GSE154891
    # retrieves info for that GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154891
    matrix_files = []
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}"
    log_message(f"fetching {url}")
    ftp = ""
    for line in execute_command(f"wget -O - {url}"):
        if m := re.search(r'.+(ftp:.*?)".+Series Matrix File', line):
            ftp = m.groups()[0]
            print(ftp)
            break
    if not ftp:
        log_message(f"No series matrix file found in {url}")
        return
    url = ftp
    ftp = ""
    log_message(f"fetching {url}")
    for line in execute_command(f"wget -O - {url}"):
        if m := re.search(r'.+(ftp:.*?)".+_series_matrix', line):
            ftp = m.groups()[0]
            print(ftp)
            break
    if not ftp:
        log_message(f"No series matrix file found in {url}")
        return
    log_message(f"fetching {ftp}")
    cmd = f"wget --quiet --no-clobber {ftp} 2> /dev/null"
    log_message(f"executing {cmd}")
    execute_command(cmd)
    file = os.path.basename(ftp)
    if os.path.exists(file):
        log_message(f"Retrieved {file}")
        matrix_files.append(file)
    else:
        log_message(f"Failed to retrieve {file} with {cmd}")
    return matrix_files

def make_row_headers_unique(data: list[list]):
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
        Sample_characteristics
        Sample_characteristics
        Sample_characteristics

    output:
        Sample_characteristics_1
        Sample_characteristics_2
        Sample_characteristics_3

    """

def bash_header(*, flags: str = "set -e"):
    return dedent(f"""
        #!/bin/bash
        {flags}
        log() {{ echo "$*" >&2 ; }}
        die() {{ echo "$*" >&2; exit 1 ; }}
    """).lstrip()


def get_ncbi_recount_data(GEO_ID: str):
    #ncbi_file = f"{GEO_ID}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    script = f"get_ncbi_recount_{GEO_ID}.sh"
    if os.path.exists(script):
        log_message(f"{script} exists for {GEO_ID}")
        return
    """
    for i in [ncbi_file, script, ncbi_file.replace(".gz", "")]:  #  spider_log_file
        if os.path.exists(i):
            log_message(f"{i} is already present for {GEO_ID}")
            return
    """
    # may not exist
    baseurl = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={GEO_ID}"
    ncbi_counts = []
    for line in execute_command(f'wget -q -O - "{baseurl}"'): # | grep NCBI.tsv.gz | grep counts_'
        if m := re.search(r'href="(.*?NCBI.tsv.gz)"', line):
            # /geo/download/?type=rnaseq_counts&amp;acc=GSE215024&amp;format=file&amp;file=GSE215024_raw_counts_GRCh38.p13_NCBI.tsv.gz
            ncbi_counts.append(m.groups()[0])
        """
        if "NCBI.tsv.gz"
        for p in line.split('"'):
            if p.endswith("NCBI.tsv.gz"):
                out = p.split("=")[-1]
                ncbi.append(f'wget -nc -O {out} "https://www.ncbi.nlm.nih.gov{p}"')
        """
    # execute_command(test_cmd, suppress_error_messages=True)
    # supply = execute_command(f'grep -i -P "length|size"  {spider_log_file}')
    if not ncbi_counts:
        log_message(f"No counts files were found for {GEO_ID} at {baseurl}")
        return
    output = ["#!/bin/bash\n"]
    for f in ncbi_counts:
        file = f.split("=")[-1]
        output.append(f"wget --no-clobber -q -O {file } 'https://www.ncbi.nlm.nih.gov{f}'")
    with open(script, "w") as tempout:
        print("\n".join(output), file=tempout)
        os.chmod(script, 0o755)
        log_message(f"source {script}")
    

def geo_parse_matrix_file(matrix_file: FileName):
    # input = name of a series matrix file
    # output = a dataframe with sample-specific info and a dataframe with common entries

    matrix_file = matrix_file.replace(".gz", "")
    gz = matrix_file + ".gz"
    if not os.path.exists(matrix_file):
        if os.path.exists(gz):
            execute_command(f"gunzip --keep {gz}")
    verify_that_paths_exist(matrix_file)
    """
    cmd = {"gz": "zcat", "txt": "cat"}.get(matrix_file.split(".")[-1], "")
    if not cmd:
        log_message(f"Unknown extension in {matrix_file} - skipping")
        return None, None
    temp = execute_command(f"{cmd} {matrix_file} | grep ^.Sample 2> /dev/null")
    if not temp:
        log_message(f"No sample info found in {matrix_file}")
        return None, None
    """
    sample_info  = []
    with open(matrix_file, "r") as tempin:
        for line in tempin:
            if line.startswith("!Sample"):
                sample_info.append(line.lstrip("!").rstrip().replace('"', "").split("\t"))

    #sample_info = [x.lstrip("!").rstrip().replace('"', "").split("\t") for x in temp]
    columns_to_drop = "Sample_relation Sample_contact_name Sample_contact_institute Sample_contact_address Sample_contact_city Sample_contact_state Sample_contact_zip/postal_code Sample_contact_country Sample_taxid_ch1 Sample_platform_id Sample_status Sample_submission_date Sample_last_update_date Sample_type Sample_data_processing Sample_growth_protocol Sample_contact_email Sample_contact_laboratory Sample_contact_department Sample_treatment_protocol Sample_extract_protocol".split()
    for row in sample_info:
        if row[0] == "Sample_channel_count":
            temp = set(row[1:])
            if len(temp) == 1 and list(temp)[0] == "1":
                columns_to_drop = [x.replace("_ch1", "") for x in columns_to_drop]
                columns_to_drop.append("Sample_channel_count")
                for i, sample in enumerate(sample_info):
                    sample_info[i][0] = sample[0].replace("_ch1", "")
            break
    sample_info = [row for row in sample_info if not (row[0] in columns_to_drop)]
    variable_columns = []
    common_columns = []
    for row in sample_info:
        values = set(row[1:])
        if len(values) == 1:
            val = list(values)[0]
            if not (val in ["None", "NONE", "0"]):
                common_columns.append(row)
                # common_columns.append(f"{row[0]}\t{val}")
        else:
            variable_columns.append(row)

    sample_info = variable_columns + common_columns
    # Sample_characteristics	cell type: iPSC
    # Sample_characteristics	condition: Control
    for i, row in enumerate(sample_info):
        if row[0] == "Sample_characteristics":
            temp = row[1].split(": ", 1)
            if len(temp) == 2:
                check = temp[0] + ": "
                if all((x.startswith(check) or x == "") for x in row[2:]):
                    sample_info[i][0] = temp[0].replace(" ", "_")
                    for j in range(1, len(row)):
                        sample_info[i][j] = sample_info[i][j].replace(check, "")

    a = make_row_headers_unique(data=sample_info)
    b = transpose_nested_list(data=a)
    return  pd.DataFrame(b[1:], columns=b[0])


def geo_find_supplementary_files(
    *, files: list[FileName], outputfile: FileName, overwrite: bool = False
):
    # input = list of series matrix files
    # output = bash script to retrieve any supplementary files, or None if None
    verify_that_paths_exist(files)
    if os.path.exists(outputfile) and overwrite is False:
        log_message(f"{outputfile} exists")
        return outputfile
    output = []
    # log_message("Ignore warnings like: Execution of wget -nc --spider ftp://... failed.")
    for matrix_file in files:
        tempoutput = []
        cmd = {"gz": "zcat", "txt": "cat"}.get(matrix_file.split(".")[-1], "")
        if cmd == "":
            log_message(f"Unknown extension in {matrix_file} - skipping")
            continue
        for temp in execute_command(
            f"{cmd} {matrix_file} | grep ^.Series_supplementary_file 2> /dev/null"
        ):
            url = temp.split('"')[1]
            if url.startswith("ftp"):
                spider_log_file = f"{os.path.basename(url)}.size"
                execute_command(
                    f"wget -nc --spider {url} 2> {spider_log_file}",
                    suppress_error_messages=True,
                )
                for g in execute_command(f'grep -i -P "length|size" {spider_log_file}'):
                    tempoutput.append(f"# {g}")
                tempoutput.append(f"wget -nc -q {url}\n")
                os.remove(spider_log_file)
            else:
                log_message(f"Malformed FTP URL in {temp}", exit_now=True)
        if tempoutput:
            output.append(f"# source: {matrix_file}")
            output += tempoutput
    if output:
        output  = [bash_header()] + output
        with open(outputfile, "w") as tempout:
            print("\n".join(output), file=tempout)
        os.chmod(outputfile, 0o755)
        log_message(f"source {outputfile} &")
        return outputfile
    else:
        return None


def geo(args):
    convert_functions = set("""
        gse-to-gsm
        gse-to-srp
        gsm-to-gse
        gsm-to-srp
        gsm-to-srr
        gsm-to-srs
        gsm-to-srx
        srp-to-gse
        srp-to-srr
        srp-to-srs
        srp-to-srx
        srr-to-gsm
        srr-to-srp
        srr-to-srs
        srr-to-srx
        srs-to-gsm
        srs-to-srx
        srx-to-gsm
        srx-to-srp
        srx-to-srr
        srx-to-srs
    """.split())

    direct_to_srp = {x.split("-to-")[0] : x for x in convert_functions if x.endswith("srp")}
    # {'srx': 'srx_to_srp', 'gsm': 'gsm_to_srp', 'gse': 'gse_to_srp', 'srr': 'srr_to_srp'}
    direct_to_gse = {x.split("-to-")[0] : x for x in convert_functions if x.endswith("gse")}
    # {'srp': 'srp_to_gse', 'gsm': 'gsm_to_gse'}

    skip = [x for x in convert_functions if x.startswith("gse")]
    convert_functions -= set(skip)
    convert_functions -= set([x for x in convert_functions if  x.split("-to-")[0] in direct_to_gse])
    # srp_to_gse  srp_to_srr  srp_to_srs  srp_to_srx  gsm_to_gse  gsm_to_srp  gsm_to_srr  gsm_to_srs  gsm_to_srx
    intermediates = set([x for x in convert_functions if  x.split("-to-")[1] in direct_to_gse])
    # srr_to_gsm     srr_to_srp    srx_to_gsm     srx_to_srp     srs_to_gsm
    indirects = set([x.split("-to-")[0] for x in intermediates])
    # {'srr', 'srx', 'srs'}
    paths_to_try = {x: [i for i in intermediates  if i.split("-to-")[0] == x] for x in indirects}
    # {'srr': ['srr_to_srp', 'srr_to_gsm'], 'srx': ['srx_to_gsm', 'srx_to_srp'], 'srs': ['srs_to_gsm']}

    #    db = SRAweb()
    #db.methods

    #obj = MyClass()
    #method_name = "my_method"
    #getattr(obj, method_name)()  # Calls obj.my_method()

    #str_to_fn = {x: eval(db.x) for  x in convert_functions}

    def call_pysradb_conversion(*, fn: str, ID: str):
       #pysradb gsm-to-gse GSM2177186
        #pysradb gsm-to-gse GSM2177186
        cmd = f"pysradb {fn} {ID}"
        log_message(f"Executing {cmd}")
        response = execute_command(f"pysradb {fn} {ID}")
        if isinstance(response, list) and len(response) > 1:
            test = response[-1].split("\t")
            if test[0] == ID:
                return test[1]
        log_message(f"empty pysradb result for {ID}")
        return ""
        
    outputfileext = args.EXT
    # outputfile = {id : f"{id}.{outputfileext}" for id in ids}
    # exit_if_files_exist(list(outputfile.values()))
    for ID in args.IDS:

        samples_output = []
        id_succession = [ID]
        src = ID.lower()[0:3]
        if src in indirects:
            for fn in paths_to_try[src]:
                if try_conversion := call_pysradb_conversion(fn=fn, ID = ID):
                    id_succession.append(try_conversion)
            if len(id_succession) > 1:
                ID = id_succession[-1]
                src = ID.lower()[0:3]
        if src in direct_to_srp:
            if new_id := call_pysradb_conversion(fn=direct_to_srp[src], ID = ID):

                id_succession.append(new_id)
                if src != "gse" and len(id_succession) > 1:
                    # don't supersede GSE to get SRP
                    ID = id_succession[-1]
                    src = ID.lower()[0:3]
        if src in direct_to_gse:
            if new_id := call_pysradb_conversion(fn=direct_to_gse[src], ID = ID):
                id_succession.append(new_id)
                ID = new_id
                src = ID.lower()[0:3]
        if src != "gse":
            log_message(f"skipping {id_succession[0]} - no GSE ID")
            continue
        get_ncbi_recount_data(GEO_ID=ID)
        outputfile = f"{ID}{outputfileext}"
        if os.path.exists(outputfile):
            log_message(f"skipping {ID} - output file {outputfile} exists")
            continue
        matrix_files = get_matrix_files_for_geo_series(geo_id=ID)

        script = f"get_supplementary_files_from_geo_{ID}.sh"
        geo_find_supplementary_files(files=matrix_files, outputfile=script)
        if problematic := [x for x in matrix_files if not ("_series_matrix" in x)]:
            temp = "\n".join(problematic)
            log_message(
                f"_series_matrix not found in all matrix files:\n{temp}", exit_now=True
            )
        subseries = [a.split("_series_matrix")[0] for a in matrix_files]
        if len(matrix_files) == 1:
            output_files = [outputfile]
        else:
            output_files = [f"{ID}.{x}.{outputfileext}" for x in subseries]
        srp_ids_found = [x for x in id_succession if x.startswith("SRP")]
        if len(srp_ids_found) == 0:
            log_message(f"No SRP ID found for {id_succession[0]}")
            SRP = ""
        else:
            if len(srp_ids_found) > 1:
                log_message(f"Multiple SRP IDs were found for {id_succession[0]}:\n" + "\n".join(srp_ids_found))
            SRP = srp_ids_found[0]
        for i, file in enumerate(matrix_files):
            metadata = geo_parse_matrix_file(file)
            metadata.insert(0, "study", id)
            metadata.insert(1, "SRP", SRP)
            metadata.insert(2, "source_", subseries[i])
            metadata = dedup_cols(data=metadata)
            metadata.to_csv(output_files[i], sep="\t", index=False)
            log_message(f"Metadata written to {output_files[i]}")
        if SRP:
            SRP = pysradb_convert_ID(fn="gse-to-srp", id=id)
            log_message(f"SRP for {id} is {SRP}")
            sra_output = f"{SRP}.sra.txt"
            if not os.path.exists(sra_output):
                #execute_command(f"{__file__} pysradb -i {SRP} -o {sra_output}")
                pysradb(SimpleNamespace(
                    EXT = "dummy", #"_pysradb.txt",
                    IDS = SRP,
                    KEEPTEMPFILES = False,
                    OUTPUTFILE = sra_output,
                    OVERWRITE = False,
                ))
            if os.path.exists(sra_output):
                sra_columns = set(
                    pd.read_csv(sra_output, sep="\t", header=0, nrows=0).columns
                )
                if matching := sra_columns & set(["experiment_alias", "library_name"]):
                    sra_column = list(matching)[0]
                    for outputfile in output_files:
                        temp = f'{outputfile.replace(".txt", "")}.{SRP}.txt'
                        cmd = f"{__file__} join -i {outputfile} {sra_output} -c Sample_geo_accession {sra_column} -o {temp}"
                        execute_command(cmd)
                        if os.path.exists(temp):
                            log_message(f"GEO + SRA => {temp}")
                        else:
                            log_message(f"Failed to create {temp} with:\n{cmd}")
                else:
                    log_message(f"No matched columns in {outputfile} and {sra_output}")

def enafastqs(args):
    # get list of ena fastqs for one or more IDs
    # called from command line
    # for a single ID, outputs list to a file or stdout
    # for multiple IDs, creates an output file for each ID
    outputfiles = {ID: f"{ID}.samples.{args.EXT}" for ID in args.IDS}
    exit_if_files_exist(list(outputfiles.values()))
    for study, outputfile in outputfiles.items():
        data = get_enafastq_list(study_ID=study, delete_temp_output=delete_temp_files)
        if os.path.exists(outputfile):
            log_message(f"{outputfile} exists - skipping {id}")
        else:
            data.to_csv(outputfile, sep="\t", index=False)
            log_message(f"{outputfile} created for {study}")
    # tempoutputfile = f"ena.{id}.temp"
    # execute_command(f"wget -O {tempoutputfile} \"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={id}&result=read_run&fields=run_accession,fastq_ftp\"")
    # https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP247494&result=read_run&fields=run_accession,fastq_ftp
    # outputfile = args.OUTPUTFILE
    """
    if len(ids) > 1:
        if outputfile:
            log_message(f"For multiple IDs, specify only {arg_def.EXT}.", exit_now=True)
        if not outputfileext:
            log_message(f"For multiple IDs, specify {arg_def.EXT}.", exit_now=True)
        outputfiles = {id : f"{id}.{outputfileext}" for id in ids}
    else:
        id = ids[0]
        outputfiles = {id : outputfile or f"{id}.{outputfileext}"}
    """


def sra_simplify_exp_desc(x):
    # GSM4503604: NSC-CB660-TERT sgTP53 + sgCDKN2A + sgPTEN +sgNF1 CLONE 1 replicate 1; Homo sapiens; RNA-Seq
    desc = str(x["experiment_title"])
    for col in set(["organism_name", "library_strategy"]) & (set(x.keys())):
        temp = "; " + str(x[col])
        desc = re.sub(temp, "", desc)
    for col in set(["library_name", "experiment_alias"]) & (set(x.keys())):
        temp = str(x[col]) + ": "
        desc = re.sub(temp, "", desc)
    desc = re.sub(" ", "_", desc)
    # log_message(f"{keep}:{desc}")
    return desc

def pysradb(args):
    test_executables("pysradb")
    ids = args.IDS
    #pysradb gsm-to-gse GSM2177186
    outputfile = args.OUTPUTFILE
    if outputfile:
        if len(ids) > 1:
            log_message(
                f"For multiple IDs, specify only --ext",
                exit_now=True,
            )
        outputfiles = {ids[0]: outputfile}
    else:
        outputfiles = {id: f"{id}{args.EXT}" for id in ids}
    overwrite = args.OVERWRITE
    if overwrite is False:
        exit_if_files_exist(list(outputfiles.values()))
    delete_temp_files = not args.KEEPTEMPFILES
    """
    if len(ids) > 1:
        if outputfile:
            log_message(f"For multiple IDs, specify only {arg_def.EXT}.", exit_now=True)
        if not ext:
            log_message(f"For multiple IDs, specify {arg_def.EXT}.", exit_now=True)
        outputfiles = {id : f"{id}.{ext}" for id in ids}
    else:
        id = ids[0]
        outputfiles = {id : outputfile or f"{id}.{ext}"}
    """
    db = SRAweb()

    for study_ID, outputfile in outputfiles.items():
        # if os.path.exists(outputfile) and delete_temp_files is False and overwrite is False:
        #    log_message(f"{outputfile} exists - skipping {id}")
        #    continue
        #data = get_sra_metadata(id=id, delete_temp_output=delete_temp_files)
        data = db.sra_metadata(study_ID, detailed=True)
        if data is None or data.shape[0] == 0:
            log_message(f"No data from SRA for {study_ID}")
            continue
        columns = list(set(["strategy", "library_strategy"]) & set(data.columns))
        if len(columns) == 0:
            log_message(f"No strategy or library_strategy for {study_ID}")
            continue
        if len(columns) > 1:
            log_message(f"Multiple strategy columns for {study_ID}")
            continue
        types = data[columns[0]].unique().tolist()
        mask = data[columns[0]] == "RNA-Seq"
        data = data[mask].copy()
        if data.empty:
            log_message(f"No RNA-seq data for {study_ID} but: " + ", ".join(types))
            continue
        for col in ["run_total_bases", "run_total_spots"]:
            data[col] = data[col].astype(int)
        data.insert(4, "read_length_temp", 0)
        data["read_length_temp"] = data["run_total_bases"] / data["run_total_spots"]
        #data["read_length_temp"] = data.apply(lambda x: x["run_total_bases"] / x["run_total_spots"], axis=1)
        data["library_layout"] = data["library_layout"].str.lower()
        columns_to_reduce = (
            "study_accession organism_name library_strategy library_layout".split()
        )
        # drop empty columns
        if empty_columns := [x for x in data.columns if all(data[x].isna()) or all(data[x] == "missing")]:
            data.drop(columns=empty_columns, inplace=True)
        columns_to_drop = """
            experiment_accession sample_accession experiment_desc library_source library_selection instrument instrument_model
            instrument_model_desc run_alias public_filename public_size public_date public_md5 public_version public_semantic_name
            public_supertype public_sratoolkit aws_url aws_free_egress aws_access_type public_url ncbi_url ncbi_free_egress 
            ncbi_access_type gcp_url gcp_free_egress gcp_access_type organism_taxid
        """.split()
        if "study_title" in data.columns:
            data["study_title"] = data["study_title"].str.replace(
                r"\s*\[RNA-Seq\]\s*", "", regex=True
            )
        if "experiment_alias" in data.columns:  # GSM
            mask = data["experiment_alias"].isna()
            tempdata = data[mask]
            if tempdata.shape[0]:
                temp = "\n".join(tempdata["run_accession"])
                log_message(f"Missing experiment_alias:\n{temp}")
            tempdata = data[~mask]
            temp = list(tempdata["experiment_alias"])
            if len(temp) != len(set(temp)):
                temp = "\n".join(temp)
                log_message(f"Non-unique experiment_alias:\n{temp}")
            if "run_alias" in data.columns:
                if list(data["run_alias"]) == list(
                    [f"{x}_r1" for x in data["experiment_alias"]]
                ):
                    columns_to_drop.append("run_alias")

        if "experiment_title" in data.columns:
            for col in set(["experiment_desc", "library_name"]) & set(data.columns):
                if list(data[col]) == list(data["experiment_title"]):
                    columns_to_drop.append(col)
            data["experiment_title"] = data.apply(
                lambda x: sra_simplify_exp_desc(x), axis=1
            )

        if columns_to_drop := set(columns_to_drop) & set(data.columns):
            data.drop(columns=list(columns_to_drop), inplace=True)
        columns_to_drop = []
        columns_to_shift = []
        # columns_to_reduce even if there are multiple values
        # drop other uninformative columns (single values)
        """
        for col in set(columns_to_reduce) & set(data.columns):
            values = list(set(data[col]))
            temp = ", ".join(values)
            output.append(f"{col}\t{temp}")
            if len(values) == 1:
                columns_to_drop.append(col)
        """
        # columns_to_reduce only if there is a single values
        for col in set(data.columns) - set(columns_to_reduce) - {"read_length_temp"}:
            values = list(set(data[col]))
            if len(values) == 1:
                # output.append(f"{col}\t{values[0]}")
                columns_to_shift.append(col)
        if ena_columns := set(
            "ena_fastq_ftp ena_fastq_ftp_1 ena_fastq_ftp_2".split()
        ) & set(data.columns):
            for col in ena_columns:
                data[col] = data[col].str.replace("era-fasp@fasp.sra.ebi.ac.uk:", "")
                # mask = data[col].str.startswith("era-fasp@fasp.sra.ebi.ac.uk:")
                # data.loc[mask, col] = data.apply(lambda x: x[col].split(":")[1], axis=1)
                # era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR214/016/SRR21488116/SRR21488116.fastq.gz
        else:
            log_message("No ENA fastqs - getting from ENA")
            ena_data = get_enafastq_list(study_ID=study_ID, delete_temp_output=delete_temp_files)
            if ena_data is None:
                log_message("No data from ENA")
            else:
                data = pd.merge(
                    data,
                    ena_data,
                    left_on=data.columns[0],
                    right_on=ena_data.columns[0],
                    how="left",
                )
                for col in ena_data.columns:
                    if all(x == "" for x in data[col]):
                        columns_to_drop.append(
                            col
                        )  # data.drop(columns = [col], inplace=True)

        if outputfile and os.path.exists(outputfile) and overwrite is False:
            # necessary to trigger deletion of temp files even if outputfile exists
            log_message(f"{outputfile} exists - skipping {study_ID}")
            continue

        data.insert(3, "read_length", 0)
        data.insert(4, "read_type", "")
        data["read_length"] = data["read_length_temp"]
        columns_to_drop.append("read_length_temp")
        # data.drop(columns = ["read_length_temp",], inplace=True)
        data.fillna("", inplace=True)
        if "ena_fastq_ftp_1" in data.columns and "ena_fastq_ftp_2" in data.columns:
            mask = (data["library_layout"] == "paired") | (
                (data["ena_fastq_ftp_1"].str.endswith(".gz"))
                & (data["ena_fastq_ftp_2"].str.endswith(".gz"))
            )
        else:
            mask = data["library_layout"] == "paired"
        data.loc[mask, "read_length"] = data.loc[mask, "read_length"] / 2
        data.loc[mask, "read_type"] = "paired"
        data.loc[~mask, "read_type"] = "single"
        data["read_length"] = data["read_length"].round(0).astype(int)
        if "paired" in set(data["read_type"]):
            data["read_length"] = data["read_length"].astype(str)
            mask = data["read_type"] == "paired"
            data.loc[mask, "read_length"] = "2x" + data.loc[mask, "read_length"]
            """
            check for paired with solo fastq
            run_accession	fastq_ftp
            SRR16106149	ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/049/SRR16106149/SRR16106149.fastq.gz
            SRR16106152	ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/052/SRR16106152/SRR16106152.fastq.gz
            SRR16106153	ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/053/SRR16106153/SRR16106153.fastq.gz
            """
            if set(["ena_fastq_ftp_1", "ena_fastq_ftp_2"]) - set(
                data.columns
            ):
                log_message("\nWarning: paired reads but single FTP file\n")
                data["sra_fastq_1"] = ""
                data["sra_fastq_2"] = ""
                mask = data["read_type"] == "paired"
                for i in [1, 2]:
                    data.loc[mask, f"sra_fastq_{i}"] = (
                        data["run_accession"] + f"_{i}.fastq"
                    )
            else:
                mask = (data["library_layout"] == "paired") & (
                    (data["ena_fastq_ftp_1"] == "") | (data["ena_fastq_ftp_2"] == "")
                )
                if mask.any():
                    log_message("\nWarning: paired reads but single FTP file\n")
                    data["sra_fastq_1"] = ""
                    data["sra_fastq_2"] = ""
                    mask = (data["library_layout"] == "paired") & (
                        (data["ena_fastq_ftp_1"] == "")
                        | (data["ena_fastq_ftp_2"] == "")
                    )
                    for i in [1, 2]:
                        data.loc[mask, f"sra_fastq_{i}"] = (
                            data["run_accession"] + f"_{i}.fastq"
                        )
        else:
            data["read_length"] = data["read_length"].astype(int)
            data["read_length"] = data["read_length"].astype(str)
        """
        for col in ["read_type", "read_length"]:
            values = ", ".join(list(set(data[col])))
            output.append(f"{col}\t{values}")
        """
        if columns_to_drop:
            data.drop(columns=list(columns_to_drop), inplace=True)
        # data = pd.DataFrame(data[1:], columns=data[0])
        if columns_to_shift:
            columns_to_shift = set(columns_to_shift)
            noshift = set(data.columns) - columns_to_shift
            new_order = [x for x in data.columns if x in noshift] + [
                x for x in data.columns if x in columns_to_shift
            ]
            data = data[new_order]  # pd.DataFrame(new_order[1:], columns=new_order[0])
        data = dedup_cols(data=data)
        data.to_csv(outputfile if outputfile else sys.stdout, sep="\t", index=False)
        """
        #output.append(outputfile else sys.stdout
        with (open(outputfile, 'w') if outputfile else sys.stdout) as tempout:
            print("\n".join(output), file=tempout)
            print("\t".join(data.columns), file=tempout)
            for i, row in data.iterrows():
                print("\t".join(map(str,row)), file=tempout)
        """
        if outputfile:
            log_message(f"{outputfile} created for {study_ID}")


def species(args):
    if args.OUTPUTFILE:
        exit_if_files_exist(args.OUTPUTFILE)
    verify_that_paths_exist(args.BAMFILES)
    species = get_species_for_bamfiles(bamfiles=args.BAMFILES)
    with open(args.OUTPUTFILE, "w") if args.OUTPUTFILE else sys.stdout as tempout:
        for bam in args.BAMFILES:
            print(f"{bam}\t{species[bam]}", file=tempout)


def get_single_species_for_bamfiles(bamfiles=None):  # , check_if_bamfiles_exist=True):
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
        temp = " ".join(species_for_bams)
        log_message(f"Invalid number of species for bam files: {temp}", exit_now=True)


def get_single_read_type_for_bamfiles(
    bamfiles=None,
):  # , check_if_bamfiles_exist=True):
    # bamfiles = list of bam files
    # reports error and stops if there's more than one read type
    # otherwise returns read type (paired or single) as a string
    assert isinstance(bamfiles, list), "bamfiles is not a list"
    # if check_if_bamfiles_exist:
    verify_that_paths_exist(bamfiles)
    # read_types = list(set([execute_command(f"{determine_read_type_sh} {bam}", exit_on_error=True) for bam in bamfiles]))
    read_types = list(
        set(
            [
                get_read_type_for_bamfile(bamfile=bam, check_if_bamfile_exists=False)
                for bam in bamfiles
            ]
        )
    )
    if len(read_types) == 1:
        return read_types[0]
    else:
        temp = " ".join(read_types)
        log_message(f"Multiple read types: {temp}", exit_now=True)

def get_genecoords(args):
    #if args.SCRIPT and args.OUTPUTFILE:
    #    log_message("Specify either --script or --outputfile but not both.", exit_now=True)
    outputfile = args.OUTPUTFILE
    outputfile = outputfile.replace("{species}", args.SPECIES)
    outputfile = outputfile.replace("{extend}", str(args.EXTEND))
    if not args.OVERWRITE:
        if outputfile:
            exit_if_files_exist(outputfile)
        if args.SCRIPT:
            exit_if_files_exist(args.SCRIPT)
    species = args.SPECIES
    coords_file = os.path.expandvars(constants.coords_source_file[species])
    verify_that_paths_exist(coords_file)
    genes_to_get = return_list_from_file_or_list(args.GENES)
    if args.SCRIPT:
        genes = " ".join(genes_to_get)
        output = f"""
            #!/bin/bash
            set -e
            ext={args.EXTEND}
            genes={genes}
            species={species}
            {__file__} get_genecoords -g $genes -s $species -x $ext -o {outputfile}
        """
        output = dedent(output)
        with open(script, "w") as tempout:
            print(output, file=tempout)
        os.chmod(script, 0o755)
        log_message(f"created {args.SCRIPT}")
    else:
        lowercase = set(x.lower() for x in genes_to_get)

        log_message(f"reading {species} gene coordinates from {coords_file}")
        data = pd.read_csv(coords_file, sep="\t", header=0)
        if " ".join(data.columns[0:4]) != "chr start stop gene_name":
            log_message(
                f"Unexpected file contents in  {coords_file}", exit_now=True
            )
        data["dummy"] = data["gene_name"].str.lower()
        mask = data["dummy"].isin(lowercase)
        data = data[mask].copy()
        found = [x for x in genes_to_get if x.lower() in set(data["dummy"])]
        if missing := "\n".join(list(set(genes_to_get) - set(found))):
            log_message(
                f"genes not found:\n{missing}"
            )  # {len(found)} genes found, {len(missing)} missing of {len(genes_to_get)} genes to get")
        data.drop(columns=["dummy"], inplace=True)
        if args.EXTEND:
            data["start"] -= args.EXTEND
            data["start"] = [
                max(1, x) for x in data["start"]
            ]  # data.apply(lambda x: max(x["start"], 1), axis=1)
            data["stop"] += args.EXTEND
            log_message("Caution: stop position may exceed chromosome length")
        data.sort_values(by=["chr", "start"]).to_csv(
            outputfile or sys.stdout, sep="\t", header=False, index=False
        )
        if outputfile:
            log_message(f"created {outputfile}")


"""
samtools fastq
  -0 FILE      write READ_OTHER to FILE
  -1 FILE      write READ1 to FILE
  -2 FILE      write READ2 to FILE
  -o FILE      write READ1 or READ2 to FILE; if a singleton file is specified with -s, only paired reads will be written to the -1 and -2 files.
  -f INT       only include reads with all  of the FLAGs in INT present [0]
  -F INT       only include reads with none of the FLAGS in INT present [0x900]
  -G INT       only EXCLUDE reads with all  of the FLAGs in INT present [0]
  -n           don't append /1 and /2 to the read name
  -N           always append /1 and /2 to the read name
  -O           output quality in the OQ tag if present
  -s FILE      write singleton READ1 or READ2 to FILE
  -t           copy RG, BC and QT tags to the FASTQ header line
  -T TAGLIST   copy arbitrary tags to the FASTQ header line, '*' for all
  -v INT       default quality score if not given in file [1]
  -c INT       compression level [0..9] to use when writing bgzf files [1]
  --i1 FILE    write first index reads to FILE
  --i2 FILE    write second index reads to FILE

Examples:
To get just the paired reads in separate files, use:
   samtools fastq -1 pair1.fq -2 pair2.fq -0 /dev/null -s /dev/null -n in.bam

To get all non-supplementary/secondary reads in a single file, redirect the output:
   samtools fastq in.bam > all_reads.fq
"""


def rnaspades(args):

    # construct commands to:
    # 1. extract unaligned reads from starting BAM (i.e., reads not mapped as proper pairs)
    # 2. extract all paired reads from extended gene-specific BAM, discarding singletons (from ends)
    # 3. run rnaspades
    # Only one fastq file is generated for each sample, but it matters for rnaspades whether the reads are single- or paired-end.
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    test_executables(constants.rnaspades.exe)

    genebamfiles = args.GENEBAMFILES or []
    genebamfiles = [x for x in genebamfiles if is_a_gene_bam(x)]
    bamfiles = args.BAMFILES or []
    bamfiles = [x for x in bamfiles if not is_a_gene_bam(x)]
    regions_file = args.REGIONS or ""

    if bamfiles and genebamfiles:
        genebam = {x: gene_bam(x) for x in bamfiles}
        if set(genebam.values()) == set(genebamfiles):
            log_message("Standard and gene BAM files match")
        else:
            log_message(
                "Problem - standard and gene BAM files specified but mismatched", exit_now=True
            )
        verify_that_paths_exist(bamfiles + genebamfiles)
    else:
        verify_that_paths_exist(bamfiles)
        genebam = {x: gene_bam(x) for x in bamfiles}
        genebamfiles = list(genebam.values())
        if any(os.path.exists(x) for x in genebamfiles):
            verify_that_paths_exist(genebamfiles)
        else:
            # need to create gene bams
            if not regions_file:
                log_message(
                    "Specify gene bams or region file to create them", exit_now=True
                )
            verify_that_paths_exist(regions_file)
            exit_if_files_exist(genebamfiles)
            for bamin, bamout in genebam.items():
                make_gene_bam(
                    inputbam=bamin, outputbam=bamout, regions=regions_file, execute =True
                )

    outputdir_by_sample = {
        x: check_dir(os.path.realpath(re.sub(r"\.bam$", "", x))) for x in bamfiles
    }

    read_type = {x: get_read_type_for_bamfile(bamfile=x) for x in bamfiles}
    fastq_files = {}
    fastq_file = {}
    for bam in bamfiles:
        prefix = re.sub(r"\.bam$", "", bam)
        if read_type[bam] == "paired":
            fastq_files[bam] = [f"{prefix}_1.fq", f"{prefix}_2.fq"]
            exit_if_files_exist(fastq_files[bam])
        else:
            fastq_file[bam] = f"{prefix}.fq"
            exit_if_files_exist(fastq_file[bam])
    cpus = args.CPUS
    mem = args.MEM

    output = []
    """
    -0 FILE      write reads designated READ_OTHER to FILE
    -1 FILE      write reads designated READ1 to FILE
    -2 FILE      write reads designated READ2 to FILE
    -o FILE      write reads designated READ1 or READ2 to FILE
                note: if a singleton file is specified with -s, only
                paired reads will be written to the -1 and -2 files.
    -f INT       only include reads with all  of the FLAGs in INT present [0]
    -F INT       only include reads with none of the FLAGS in INT present [0x900]
    -G INT       only EXCLUDE reads with all  of the FLAGs in INT present [0]
    -n           don't append /1 and /2 to the read name
    -N           always append /1 and /2 to the read name
    -s FILE      write singleton reads designated READ1 or READ2 to FILE
    """
    for bam in bamfiles:
        output.append(f"\n# {bam}")
        outputdir = outputdir_by_sample[bam]
        log = f"{outputdir}.log"
        if read_type[bam] == "paired":
            # library consists of paired reads - extract reads not aligned by STAR as proper pairs
            fq1, fq2 = fastq_files[bam]
            temp1 = f"{genebam[bam]}.temp.1"
            temp2 = f"{genebam[bam]}.temp.2"
            not_aligned_as_pair = "-F 0x2"  # not aligned as proper pair, though STAR results appear to contain only proper pairs or unaligned
            output.append(dedent(f"""
                fq1={fq1}
                fq2={fq2}
                temp1={temp1}
                temp2={temp2}
                if [ -e $fq1 ] && [ -e $fq2 ] && [ -e $temp1 ] && [ -e $temp2 ] ; then
                    log "all input files exist"
                else
                    samtools fastq {not_aligned_as_pair} --threads {cpus} -1 $fq1 -2 $fq2 -N -s /dev/null -0 /dev/null {bam}
                    lines1=$(wc -l $fq1 | cut -f 1 -d " ")
                    lines2=$(wc -l $fq2 | cut -f 1 -d " ")
                    log "$fq1: $lines1 lines"
                    log "$fq2: $lines2 lines"
                    if [ $lines1 != $lines2 ] ; then die "mismatch between $fq1 and $fq2" ; fi
                    # 
                    samtools sort -n -O SAM {genebam[bam]} | samtools fastq -1 $temp1 -2 $temp2 -s /dev/null -N -0 /dev/null
                    lines1=$(wc -l $temp1 | cut -f 1 -d " ")
                    lines2=$(wc -l $temp2 | cut -f 1 -d " ")
                    log "$temp1: $lines1 lines"
                    log "$temp2: $lines2 lines"
                    if [ $lines1 != $lines2 ] ; then die "mismatch between $temp1 and $temp2" ; fi
                    cat $temp1 >> $fq1 && cat $temp2 >> $fq2
                fi
                # spades
                {constants.rnaspades.exe} -1 $fq1 -2 $fq2 -t {cpus} -m {mem} -o {outputdir} &>> {log}
            """).lstrip())
            """
            output.append(
                f"if [ -e {fq1} ] && [ -e {fq2} ] && [ -e {temp1} ] && [ -e {temp2} ] ; then"
            )
            # output.append(f"if [[ -e {fq1} ]] ; then")
            output.append(f'  echo "# input files exist "')
            output.append("else")

            output.append(
                f"  samtools fastq {not_aligned_as_pair} --threads {cpus} -1 {fq1} -2 {fq2} -N -s /dev/null -0 /dev/null {bam}"
            )
            output.append(f'  lines1=$(wc -l {fq1} | cut -f 1 -d " ")')
            output.append(f'  lines2=$(wc -l {fq2} | cut -f 1 -d " ")')
            output.append(f'  echo "# {fq1}: $lines1 lines"')
            output.append(f'  echo "# {fq2}: $lines2 lines"')
            output.append(
                f'  if [ $lines1 != $lines2 ] ; then echo "mismatch between {fq1} and {fq2}" ; exit ; fi'
            )
            output.append(
                f"  samtools sort -n -O SAM {genebam[bam]} | samtools fastq -1 {temp1} -2 {temp2} -s /dev/null -N -0 /dev/null"
            )
            output.append(f'  lines1=$(wc -l {temp1} | cut -f 1 -d " ")')
            output.append(f'  lines2=$(wc -l {temp2} | cut -f 1 -d " ")')
            output.append(f'  echo "# {temp1}: $lines1 lines"')
            output.append(f'  echo "# {temp2}: $lines2 lines"')
            output.append(
                f'  if [ $lines1 != $lines2 ] ; then echo "mismatch between {temp1} and {temp2}" ; exit ; fi'
            )
            output.append(f"  cat {temp1} >> {fq1} && cat {temp2} >> {fq2}")
            output.append("fi")
            # spades
            # == Error ==  file is empty: ${HOME}/star/kdm6a.chan/RT4_KO_2.genes.bam.temp.singletons.fq (single reads, library number: 1, library type: paired-end)
            output.append(
                f"{constants.rnaspades.exe} -1 {fq1} -2 {fq2} -t {cpus} -m {mem} -o {outputdir} &>> {log}"
            )
            """
            delete = [fq1, fq2, temp1, temp2]
        else:
            # library consists of single reads - extract reads not aligned by STAR
            unaligned = "-f 0x4"  # 0x4     4  UNMAP          segment unmapped
            delete = fq = fastq_file[bam]
            temp = os.path.join(outputdir, os.path.basename(f"{genebam[bam]}.temp"))
            outputfastq = os.path.join(outputdir, os.path.basename(fq))
            output.append(dedent(f"""
                if [[ -e {fq} ]] ; then
                    log "{fq} exists"
                else
                    samtools fastq {unaligned} --threads {cpus} {bam} > {fq} 2> {log}
                    samtools fastq {genebam[bam]} >> {fq} 2>> {log}
                fi
                {constants.rnaspades.exe} -s {fq} -t {cpus} -m {mem} -o {outputdir} &>> {log}
            """).lstrip())

            """
            output.append(f"if [[ -e {fq} ]] ; then")
            output.append(f'  log "# {fq} exists"')
            output.append("else")
            output.append(
                f"  samtools fastq {unaligned} --threads {cpus} {bam} > {fq} 2> {log}"
            )
            output.append(f"  samtools fastq {genebam[bam]} >> {fq} 2>> {log}")
            output.append(f"fi")
            output.append(
                f"{constants.rnaspades.exe} -s {fq} -t {cpus} -m {mem} -o {outputdir} &>> {log}"
            )
            """

            delete = [fq]  #
        delete = " ".join(delete)
        output.append(dedent(f"""
            # cleanup
            rm -fr {delete} {outputdir}/transcripts.paths {outputdir}/misc \\
                {outputdir}/K* {outputdir}/*fastg {outputdir}/*gfa {outputdir}/before_rr.fasta \\
                {outputdir}/pipeline_state {outputdir}/tmp
            nice gzip -q {outputdir}/*fasta &
        """))
        if args.ZIPBAMS:
            output.append(dedent(f"""
                # zip gene BAM files
                zip -qymT -u {constants.gene_bams_zipfile} {genebam[bam]} {genebam[bam]}.bai 2> /dev/null"

            """))

    outputdir = os.path.dirname(outputdir)
    output.append(dedent(f"""
        # find {outputdir} -type l -delete 2> /dev/null")
        find . -name *.log 2> /dev/null | zip -mT logs.zip 2> /dev/null
        find . -name params.txt 2> /dev/null | zip -mT -@ spades.config.zip 2> /dev/null
        find . -name dataset.info 2> /dev/null | zip -mT -@ spades.config.zip 2> /dev/null
        find . -name *.yaml 2> /dev/null | zip -mT -@ spades.config.zip 2> /dev/null
    """))

    # output = ["#!/bin/bash\n"] + output # avoid set -e
    cmd = f"{outputscript} &> {outputscript}.log &" if outputscript else ""
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )


def handle_commands_and_output(*,
    commands: str | list[str], outputfile: FileName|None = None, single_line_command: str|None =  None, execute: bool = False
):
    # commands = list of strings
    # single_line_command = individual command to execute script, only valid if outputfile is specified
    # exec: execute single_line_command
    assert (
        single_line_command is None or outputfile is not None
    ), "single_line_command requires outputfile"
    assert (
        execute is False or single_line_command is not None
    ), "execute is True but single_line_command is None"
    if outputfile:
        if isinstance(commands, str):
            commands = [commands]
        if not (commands[0].startswith("#!/bin/bash")):
            commands = [bash_header()] + commands
            log_message(f"Bash header missing in call from {_caller()}")
        with open(outputfile, "w") as tempout:
            print("\n".join(commands), file=tempout)
        # execute_command(f"chmod +x {outputfile}")
        os.chmod(outputfile, 0o755)
        log_message(f"# created {outputfile}")
        if single_line_command:
            if execute:
                execute_command(single_line_command)
            else:
                print(single_line_command)
    else:
        if execute:
            for cmd in commands:
                log_message(f"# executing {cmd}")
                execute_command(cmd)
        else:
            print("\n".join(commands))


def make_gene_bam(*, inputbam: FileName, outputbam: FileName, regions: FileName, execute: bool = False):
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


def genebams(args):
    # test_executables(exes="samtools"])
    bamfiles = [x for x in args.BAMFILES if is_a_gene_bam(x) == False]
    verify_that_paths_exist(bamfiles + [f"{bam}.bai" for bam in bamfiles])
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    regions_file = args.REGIONS
    verify_that_paths_exist(regions_file)
    outputbam = {x: gene_bam(x) for x in bamfiles}
    exit_if_files_exist(list(outputbam.values()))
    output = []
    if args.EXEC:
        for bamin, bamout in outputbam.items():
            make_gene_bam(
                inputbam=bamin, outputbam=bamout, regions=regions_file, execute =True
            )
    else:
        for bamin, bamout in outputbam.items():
            output.append(
                make_gene_bam(
                    inputbam=bamin, outputbam=bamout, regions=regions_file, execute =False
                )
            )

    if args.ZIPBAMS:
        zipcommands = [
            f"zip -qymT -u {constants.gene_bams_zipfile} {bam} {bam}.bai"
            for bam in outputbam.keys()
        ]
        if exec:
            for x in zipcommands:
                execute_command(command=x, exit_on_error=True, log_command=True)
        else:
            output += zipcommands
    if output:
        cmd = f"source {outputscript} &> {outputscript}.log" if outputscript else ""
        handle_commands_and_output(
            commands=output,
            outputfile=outputscript,
            single_line_command=cmd,
            execute =args.EXEC,
        )
    # else:
    #    log_message("Error - nothing to do", exit_now=True)


def junctions(args):
    from tempfile import mkstemp

    # needs coords and either gene bams or non-gene BAMs, but not both
    test_executables("regtools")
    outputfileext = args.EXT
    bamfiles = args.BAMFILES
    verify_that_paths_exist(bamfiles)
    #species = args.SPECIES or get_single_species_for_bamfiles(bamfiles=bamfiles)
    #cpus = args.CPUS
    outputfile = {bamfile: f"{bamfile}.{outputfileext}" for bamfile in bamfiles}
    # outputfileprefix = {x: re.sub('\.bam$', f".realigned", x) for x in genebamfiles}
    exit_if_files_exist(outputfile.values())
    mincount = args.MINCOUNT
    if len(set(outputfile.values())) != len(set(outputfile.keys())):
        log_message("Redundant output file(s)", exit_now=True)
    regions = args.REGIONS
    if regions:
        verify_that_paths_exist(regions)
        regiondata = pd.read_csv(regions, sep="\t", header=None)
        regiondata = regiondata[regiondata.columns[0:3]]
        regiondata.columns = ["chr", "start", "stop"]
    tempoutputfile = mkstemp(dir=os.getcwd())
    for bamfile in bamfiles:
        if regions:
            for i, row in regiondata.iterrows():
                cmd = f'regtools junctions extract -s RF -r {row["chr"]}:{row["start"]}-{row["stop"]} {bamfile} >> tempoutputfile'
                execute_command(cmd)
        else:
            cmd = f"regtools junctions extract -s RF -o {tempoutputfile} {bamfile}"
            execute_command(cmd)
        data = pd.read_csv(
            tempoutputfile, sep="\t", header=None, names=constants.bed12_columns
        )
        data["strand"] = "?"
        data["name"] = ""
        data["_score"] = data["score"]
        data["score"] = 0
        data = data.groupby(constants.bed12_columns, as_index=False).sum()
        data["score"] = data["_score"]
        if mincount:
            mask = data["score"] >= mincount
            data = data[mask].copy()
        data["name"] = [f"JUNC{n+1:08d}" for n in range(data.shape[0])]
        data.drop(columns=["_score"]).to_csv(outputfile[bamfile], sep="\t", index=False)
    os.remove(tempoutputfile)


"""
		-a INT	Minimum anchor length. Junctions which satisfy a minimum anchor length on both sides are reported. [8]
		-m INT	Minimum intron length. [70]
		-M INT	Maximum intron length. [500000]
		-o FILE	The file to write output to. [STDOUT]
"""


def abra2(args):
    # needs coords and either gene bams or non-gene BAMs, but not both
    test_executables("abra2")
    genebamfiles = args.GENEBAMFILES or []
    genebamfiles = [x for x in genebamfiles if is_a_gene_bam(x)]
    bamfiles = args.BAMFILES or []
    bamfiles = [x for x in bamfiles if not is_a_gene_bam(x)]

    outputscript = args.SCRIPT or ""
    if outputscript:
        exit_if_files_exist(outputscript)
    regions_file = args.REGIONS
    verify_that_paths_exist(regions_file)
    cpus = args.CPUS - 1
    if bamfiles:
        genebam_for_bam = {x: gene_bam(x) for x in bamfiles}
        if genebamfiles:
            if set(genebam_for_bam.values()) == set(genebamfiles):
                log_message("Standard BAMs ignored since gene bams exist")
                bamfiles = []
            else:
                log_message(
                    "Problem - bams & gene bams specified but mismatched", exit_now=True
                )
    if bamfiles:
        # make gene bams
        verify_that_paths_exist(bamfiles)
        genebamfiles = list(genebam_for_bam.values())
        exit_if_files_exist(genebamfiles)
        for bamin, bamout in genebam_for_bam.items():
            make_gene_bam(
                inputbam=bamin, outputbam=bamout, regions=regions_file, execute =True
            )
    if not genebamfiles:
        log_message("No input BAM or gene BAM file(s)", exit_now=True)
    verify_that_paths_exist(genebamfiles)
    outputbam = {x: realigned_gene_bam(x) for x in genebamfiles}
    exit_if_files_exist(list(outputbam.values()))
    species = args.SPECIES or get_single_species_for_bamfiles(bamfiles=genebamfiles)
    read_type = args.READS or get_single_read_type_for_bamfiles(bamfiles=genebamfiles)
    command = f"abra2 --nosort --threads {cpus} --undup"
    if read_type == "single":
        command += " --single"
    elif read_type != "paired":
        log_message(f"Invalid read type {read_type}", exit_now=True)
    regiondata = pd.read_csv(regions_file, sep="\t", header=None, dtype=object)
    chr_column = regiondata.columns[0]
    chrs = sorted(set(regiondata[chr_column]))

    ref_dir = os.path.join(constants.ref_dir, f"by.chr.{species}")
    ref_file_prefix_by_chr = {chr: os.path.join(ref_dir, chr) for chr in chrs}
    files_not_found = []
    output = []
    output.append(bash_header())
    for chr in chrs:
        for ext in [".gtf", ".fa"]:
            ref_file = ref_file_prefix_by_chr[chr] + ext
            if os.path.exists(ref_file):
                pass
            elif os.path.exists(f"{ref_file}.gz"):
                output.append(f"gunzip -q {ref_file}.gz")
            else:
                files_not_found.append(ref_file)
    if errors := "\n".join(files_not_found):
        log_message(f"ref files not found:\n{errors}", exit_now=True)

    outputfileprefix = {x: re.sub(r"\.bam$", f".realigned", x) for x in genebamfiles}
    temp_region_file_by_chr = {chr: f"{regions_file}.temp.{chr}" for chr in chrs}
    if len(set(temp_region_file_by_chr.values())) != len(
        set(temp_region_file_by_chr.keys())
    ):
        log_message(f"length mismatch {str(temp_region_file_by_chr)}", exit_now=True)

    inputlist = ",".join(genebamfiles)
    for chr, tempdata in regiondata.groupby(chr_column):
        if chr.startswith("chr"):
            label = chr
        else:
            label = f"chr{chr}"
        output.append(f"\n# {label}")
        tempdata.to_csv(
            temp_region_file_by_chr[chr], sep="\t", header=False, index=False
        )
        chrgtf = ref_file_prefix_by_chr[chr] + ".gtf"
        chrfa = ref_file_prefix_by_chr[chr] + ".fa"
        reg = temp_region_file_by_chr[chr]
        outputlist = ",".join(
            [f"{outputfileprefix[bam]}.temp.{label}.bam" for bam in genebamfiles]
        )
        log = f"abra2.{label}.log"
        output.append(
            f"{command} --gtf {chrgtf} --ref {chrfa} --targets {reg} --in {inputlist} --out {outputlist} &>> {log}"
        )

    for inputbam, prefix in outputfileprefix.items():
        output.append(f"\n# {inputbam}")
        newbam = f"{prefix}.bam"
        output.append(
            f"samtools cat -o - {prefix}.temp.*.bam | samtools sort -o {newbam} -"
        )
        output.append(f"rm {prefix}.temp.*.bam")
        output.append(f"samtools index {newbam}")
        if args.ZIPBAMS:
            output.append(
                f"zip -qymT -u {constants.gene_bams_zipfile} {inputbam} {inputbam}.bai 2> /dev/null"
            )
            output.append(
                f"zip -qymT -u {constants.gene_bams_realigned_zipfile} {newbam} {newbam}.bai 2> /dev/null"
            )
    temp = " ".join(list(temp_region_file_by_chr.values()))
    output.append(
        f"zip -qymT abra2.zip {outputscript} {temp} {regions_file} abra2.*.log 2> /dev/null"
    )
    output.append(f"# find . -type l -delete")
    output.append(f"# find . -maxdepth 2 -empty -delete")

    cmd = f"nohup {outputscript} &>> {outputscript}.log &" if outputscript else ""
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )

def return_list_from_file_or_list(genes: list[str]):
    if len(genes) == 1:
        file = os.path.expandvars(genes[0])
        if os.path.exists(file):
            return list(pd.read_csv(file, sep="\t", header=None, dtype=object)[0])
    return genes


def return_dict_from_file_or_list(args=None, arglabel=None):
    assert arglabel in args.keys(), f"{arglabel} not found in args.keys()"
    argvalues = args[arglabel]
    assert isinstance(argvalues, list), f"args[{arglabel}] is not a list"
    if len(argvalues) == 1 and os.path.exists(argvalues[0]):
        file = argvalues[0]
        # log_message(f"Obtaining {arglabel} from {file}")
        temp = pd.read_csv(file, sep="\t", header=None, dtype=object)
        if len(temp.columns) < 2:
            log_message(f"Fewer than two columns found in {file}", exit_now=True)
        return (
            temp[temp.columns[0:2]].set_index(temp.columns[0]).T.to_dict("records")[0]
        )
    # log_message(f"Obtaining {arglabel} from command line")
    # expecting from:to pairs
    rename = {}
    for pair in argvalues:
        temp = pair.split(":")
        if len(temp) != 2:
            log_message(f"Malformed from:to pair {pair}", exit_now=True)
        rename[temp[0]] = temp[1]
    return rename


def find_gtf(species=None):
    gtf = constants.gtf.get(species, "")
    if not gtf:
        log_message(f"Undefined default gtf for {species}", exit_now=True)
    if os.path.exists(gtf):
        return gtf
    if os.path.exists(f"{gtf}.gz"):
        log_message(f"zcat {gtf}.gz > {gtf}\nor\ngunzip {gtf}", exit_now=True)
    log_message(f"{gtf} not found for {species}", exit_now=True)


def featureCounts(args):
    # expects all BAM files to be from the same species
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    bamfiles = args.BAMFILES
    verify_that_paths_exist(bamfiles)

    sort_by_read_ID = args.SORT_BY_READ_ID == "yes"
    if sort_by_read_ID:
        unsortedbam = {bam: re.sub(r"\.bam$", ".unsorted.bam", bam) for bam in bamfiles}
        # already_sorted = [bam for bam in bamfiles if os.path.exists(unsortedbam[bam])]
        outputfile = {
            bam: re.sub(r"\.bam$", ".unsorted.counts.txt", bam) for bam in bamfiles
        }
    else:
        outputfile = {bam: re.sub(r"\.bam$", ".counts.txt", bam) for bam in bamfiles}
    exit_if_files_exist(outputfile.values())
    species = args.SPECIES or get_single_species_for_bamfiles(bamfiles=bamfiles)

    # constants.sortcpus = floor(args.CPUS / 2)

    
    # log_message(f"constants.sortcpus: {constants.sortcpus}\nconstants.sortmem: {constants.sortmem}")
    if (args.SORT_CPUS + 1) * args.SORTMEM > system_memory():
        log_message(
            f"{args.SORT_CPUS +1} CPUs x {args.SORTMEM}M > {system_memory()}", exit_now=True
        )
    ref_file = find_gtf(species=species)

    output = []
    if outputscript:
        output.append(bash_header())
    # if unzip_commands:
    #    output.append("\n".join(unzip_commands))
    for bam in bamfiles:
        gtf = constants.gtf[species_per_bam[bam]]
        out = outputfile[bam]
        if sort_by_read_ID:
            # works keep
            tempbam = unsortedbam[bam]
            output.append(dedent(f"""
                samtools view -F 0x4 --bam {bam} 2> {out}.unsort.1.err | samtools sort -n -@ {args.SORT_CPUS} -m {args.SORTMEM}M -l 9 -o {tempbam} 2> {out}.unsort.2.err
                featureCounts {constants.featurecounts.options} -a {gtf} -o {out} -T {args.CPUS} {tempbam} 2> {out}.log
            """)).lstrip()
            if outputscript:
                output.append(f"rm {tempbam}")

            # avoid writing temp BAM file
            # print(f"", file=tempout)
            # print(f"samtools view -F 0x4 --bam {bam} 2> {out}.view.err | samtools sort -n -@ {constants.sortcpus} -m {constants.sortmem} -O sam 2> {out}.sort.err | featureCounts -a {gtf} -o {out} {opts} -T {constants.sortcpus} 2> {out}.log", file=tempout)
        else:
            # this would fail since opts is unspecified
            output.append(f"featureCounts -a {gtf} -o {out} {opts} {bam} 2>> {out}.log")
    with open(outputscript, "w") if outputscript else sys.stdout as tempout:
        print("\n".join(output), file=tempout)

    cmd = f"nohup {outputscript} &>> {outputscript}.log &" if outputscript else ""
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )


def unsort(args):
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    bamfiles = args.BAMFILES
    verify_that_paths_exist(bamfiles)
    outputfile = {bam: re.sub(r"\.bam$", ".unsorted.bam", bam) for bam in bamfiles}
    exit_if_files_exist(outputfile.values())

    sortcpus = args.SORT_CPUS
    sortmem = args.SORTMEM
    log_message(
        f"sortcpus: {sortcpus}\nsortmem: {sortmem}"
    )
    if (sortcpus + 1) * sortmem > system_memory():
        log_message(
            f"{sortcpus +1} CPUs x {sortmem}M > {system_memory()}",
            exit_now=True,
        )
    output = []
    output.append(bash_header())
    for bam in bamfiles:
        out = outputfile[bam]
        # print(f"samtools view -F 0x4 --bam -@ {cpus} {bam} 2> {out}.1.err | samtools sort -n -@ {sortcpus} -m {sortmem}M -l 9 -o {out} 2> {out}.2.err", file=tempout)
        output.append(
            f"samtools view -F 0x4 --bam {bam} 2> {out}.1.err | samtools sort -n -@ {sortcpus} -m {sortmem}M -l 9 -o {out} 2> {out}.2.err"
        )

    cmd = f"nohup {outputscript} &>> {outputscript}.log &" if outputscript else ""
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )


def unique_gene_metadata(x):
    return ";".join(list(set(str(x).split(";"))))

def gene_metadata(args):
    log_message("output gene metadata")
    # inputfile = args.INPUTFILES[0]
    verify_that_paths_exist(args.INPUTFILE)
    outputfile = args.OUTPUTFILE
    if outputfile:
        if (b := os.path.basename(outputfile)) != outputfile:
            log_message(f"outputfile {outputfile} != basename {b}", exit_now=True)
        outputdir = args.OUTPUTDIR.rstrip("/")
        outputfile = os.path.join(outputdir, outputfile)
        exit_if_files_exist(outputfile)
    # output gene metadata from featurecounts
    verify_that_paths_exist(input)
    data = pd.read_csv(input, sep="\t", header=0, dtype=object, comment="#", nrows=1)
    columns_to_keep = list(
        set(constants.featurecounts.constant_columns) & set(data.columns)
    )
    # keep = list(data.columns)
    # keep.pop()
    data = pd.read_csv(
        input, sep="\t", header=0, dtype=object, comment="#", usecols=columns_to_keep
    )
    # data.applymap(lambda x: ";".join(list(set(x.split(";")))))
    for col in data.columns[1:]:
        data[col] = data[col].apply(lambda x: unique_gene_metadata(x))
    if output:
        data.to_csv(output, sep="\t", index=False)
        log_message(f"Gene metadata written to {output}")
    else:
        return data



def genomeref(args):
    if args.SCRIPT:
        if "{species}" in args.SCRIPT:
            args.SCRIPT = args.SCRIPT.replace("{species}", args.SPECIES)
        if not args.OVERWRITE:
            exit_if_files_exist(args.SCRIPT)
    dna = constants.ref.dna[args.SPECIES]
    rna = constants.ref.rna[args.SPECIES]
    output = []
    output.append(bash_header())
    output.append(dedent(f"""
        destdir={args.OUTPUTDIR}
        dna={dna}
        rna={rna}
        wget --no-clobber -P $destdir $dna
        wget --no-clobber -P $destdir $rna
    """))
    if args.SCRIPT:
        with open(args.SCRIPT, "w") as tempout:
            print("\n".join(output), file=tempout)
        os.chmod(args.SCRIPT, 0o755)
        log_message(f"created {args.SCRIPT}")
    else:
        print("\n".join(output))


"""
input if source is aspera (tabs only, header optional and arbitrary)
    sample	ena_fastq_ftp_1	                                            ena_fastq_ftp_2
    WT_1	vol1/fastq/SRR116/090/SRR11624290/SRR11624290_1.fastq.gz	vol1/fastq/SRR116/090/SRR11624290/SRR11624290_2.fastq.gz
    WT_2	vol1/fastq/SRR116/091/SRR11624291/SRR11624291_1.fastq.gz	vol1/fastq/SRR116/091/SRR11624291/SRR11624291_2.fastq.gz
or:
    sample	ENA_aspera
    WT_1	vol1/fastq/SRR296/005/SRR2966825/SRR2966825.fastq.gz
    WT_2	vol1/fastq/SRR296/006/SRR2966826/SRR2966826.fastq.gz

input if  source = sra (tabs only, header optional and arbitrary)
    sample	    fastq_1	            fastq_2
    WT_1	    SRR11624290_1.fastq	SRR11624290_2.fastq
    KO_TP53_1	SRR11624292_1.fastq	SRR11624292_2.fastq
or:
    sample	        fastq
    WT_1	        SRR2966825.fastq
    WT_2	        SRR2966826.fastq
"""


# run before featurecounts if starting from coordinate-sorted BAM
def bash_samtools_sort_by_read_ID(sortcpus=None, sortmem=None):
    bash_function_name = "sort_bam_by_coord"
    code = f"""
        {bash_function_name}(){{
            inputbam=$1
            outputbam=$2
            err1=$outputbam.unsort.1.log
            err2=$outputbam.unsort.2.log
            samtools view -F 0x4 --bam $inputbam 2> $err1 | samtools sort -n -@ {sortcpus} -m {sortmem}M -l 8 -o $outputbam 2> $err2
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


# run after STAR and featurecounts
def bash_samtools_sort_by_position(sortcpus=None, sortmem=None):
    bash_function_name = "sort_bam_by_coord"
    code = f"""
        {bash_function_name}(){{
            inputbam=$1
            outputbam=$2
            #log=$outputbam.sort.log
            samtools sort -@ {sortcpus} -m {sortmem}M -l 8 -o $outputbam $inputbam #2> $log
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def bash_featurecounts(*, gtf: FileName, cpus: int, options: str):
    opts = f"{options} -T {cpus}"
    bash_function_name = "run_featurecounts"
    code = f"""
        {bash_function_name}(){{
                bamfile=$1
                outputfile=$2
                gtf={gtf}
                options="{opts}"
                featureCounts $options -a $gtf -o $outputfile  $bamfile 2> $outputfile.log
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def bash_aspera(speed=0):
    bash_function_name = "fastqs_aspera"
    code = f"""
        {bash_function_name}(){{
            file_list=$1
            bin=${HOME}/.aspera/connect/bin/ascp
            ssh=${HOME}/.aspera/connect/ls asperaweb_id_dsa.openssh
            opts="--mode recv --user era-fasp --host fasp.sra.ebi.ac.uk --overwrite=never -QT -P 33001"
            speed="-l {speed}m"    
            $bin $speed $opts -i $ssh --file-list=$file_list . 2> $file_list.log
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def bash_fasterq(*, cpus: int, read_type: define_constants().read_types):
    bash_function_name = "fastqs_sra"
    opts = {"single": "", "paired": "--split-files"}.get(read_type)
    code = f"""
        {bash_function_name}(){{
            opts="{opts} --threads {cpus} --bufsize 500MB --curcache 500MB --mem 2500MB"
            for SRRid in "$@" ; do
            if [[ -e $SRRid ]] ; then
                log "skipping $SRRid"
            else
                fasterq-dump $opts $SRRid 2> $SRRid.log
            fi
            done
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def bash_star(*, cpus: int = 1, index=None, readFilesCommand: ["cat", "zcat"]):
    assert readFilesCommand in ["cat", "zcat"], f"Invalid readFilesCommand {readFilesCommand}"
    bash_function_name = "run_star"
    code = f"""
        {bash_function_name}(){{
            fastqs=$1
            prefix=$2
            otheropts=$3
            out=$prefix.star.out
            err=$prefix.star.err
            sampleopts="--outFileNamePrefix $prefix."
            inputopts="--readFilesCommand {readFilesCommand}"
            index={index}
            outputopts="--outSAMunmapped Within --outBAMcompression 8 --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 0 --outSAMtype BAM Unsorted --runDirPerm All_RWX"
            STAR --genomeLoad LoadAndKeep --runMode alignReads --genomeDir $index --runThreadN {cpus} \\
                --readFilesIn $fastqs 
                $inputopts \\
                $outputopts \\
                $otheropts \\
                $sampleopts  > $out 2> $err
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def parse_fastq_manifest(
    inputfile=None, fastq_source=None, outputdir=None, overwrite=False
):
    # read fastq manifest, common checks
    assert fastq_source in ["ENA", "SRA"], f"invalid fastq source {fastq_source}"
    verify_that_paths_exist(inputfile)
    data = pd.read_csv(inputfile, sep="\t", header=None, dtype=object, nrows=1)
    if len(data.columns) < 2:
        log_message(f"Insufficient columns in {inputfile}", exit_now=True)
    if len(data.columns) > 3:
        log_message(f"Too many columns found in {inputfile}", exit_now=True)
    if ".fastq" in data[data.columns[1]][0] or data[data.columns[1]][0].startswith(
        "SRR"
    ):
        header = None
        log_message(f"No header row found in {inputfile}")
    else:
        header = 0
        for i, row in data[0:1].iterrows():
            temp = "\t".join(row)
            log_message(f"Header row found in {inputfile}:\n{temp}")
    data = pd.read_csv(inputfile, sep="\t", header=_header, dtype=object)
    errors = False
    label = data.columns[0]
    label = data.columns[0]
    samples = data[label].unique().tolist()

    if problematic := [sample for sample in samples if " " in sample]:
        temp = "\n".join(problematic)
        log_message(f"Spaces in sample IDs:\n{temp}")
        errors = True
    if problematic := [sample for sample in samples if ".fastq" in sample]:
        temp = "\n".join(problematic)
        log_message(
            f".fastq found in first column, which should be a sample label:\n{temp}"
        )
        errors = True
    fastq_columns = list(data.columns[1:])
    if len(fastq_columns) == 1:
        read_type = "single"
    elif len(fastq_columns) == 2:
        read_type = "paired"
    else:
        log_message("Invalid number of columns", exit_now=True)
    log_message(f"{read_type} reads")
    if fastq_source == "ENA":
        for col in fastq_columns:
            data[col] = data[col].str.replace(
                "era-fasp@fasp.sra.ebi.ac.uk:", "", regex=False
            )
        if problematic := [
            x
            for x in data[col]
            for col in fastq_columns
            if not re.match(r"^vol1/.+fastq.gz$", x)
        ]:
            temp = "\n".join(problematic)
            log_message(f"fastq files not matching ^vol1/.+fastq.gz$:\n{temp}")
            errors = True
        if read_type == "single":
            if problematic := [
                x for x in data[fastq_columns[0]] if re.search(r"_1.fastq|_2.fastq", x)
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid fastq files for {read_type} reads:\n{temp}")
                errors = True
        else:
            if problematic := [
                x for x in data[fastq_columns[0]] if not "_1.fastq" in x
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid 1 fastq files for {read_type} reads:\n{temp}")
                errors = True
            if problematic := [
                x for x in data[fastq_columns[1]] if not "_2.fastq" in x
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid 2 fastq files for {read_type} reads:\n{temp}")
                errors = True
    elif fastq_source == "SRA":
        if read_type == "single":
            if problematic := [
                x for x in data[fastq_columns[0]] if not re.match(r"^SRR\d+.fastq$", x)
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid fastq files for {read_type} reads:\n{temp}")
                errors = True
        else:
            if problematic := [
                x
                for x in data[fastq_columns[0]]
                if not re.match(r"^SRR\d+_1.fastq$", x)
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid 1 fastq files for {read_type} reads:\n{temp}")
                errors = True
            if problematic := [
                x
                for x in data[fastq_columns[1]]
                if not re.match(r"^SRR\d+_2.fastq$", x)
            ]:
                temp = "\n".join(problematic)
                log_message(f"invalid 2 fastq files for {read_type} reads:\n{temp}")
                errors = True
    if errors:
        log_message("Errors found in fastq file manifest", exit_now=True)

    star_fastq_files = {x: "" for x in samples}
    fastq_fetch_params = {x: "" for x in samples}

    if fastq_source == "ENA":
        fastq_file_list = {
            sample: os.path.join(outputdir, f"fastq_files_aspera.{sample}")
            for sample in samples
        }
        # output fastq lists
        if not overwrite:
            exit_if_files_exist(list(fastq_file_list.values()))
        fastq_fetch_params = fastq_file_list

    for sample in samples:
        mask = data[label] == sample
        samplefastqs = data[mask]
        if fastq_source == "ENA":
            with open(fastq_file_list[sample], "w") as tempout:
                for col in fastq_columns:
                    print("\n".join(samplefastqs[col]), file=tempout)
                    temp = map(os.path.basename, samplefastqs[col])
                    star_fastq_files[sample] += ",".join(temp) + " "
        else:
            # SRA
            for col in fastq_columns:
                star_fastq_files[sample] += ",".join(samplefastqs[col]) + " "
            srrIDs = [
                x.split("_")[0].split(".")[0] for x in samplefastqs[fastq_columns[0]]
            ]
            fastq_fetch_params[sample] = " ".join(srrIDs)
        star_fastq_files[sample] = star_fastq_files[sample].rstrip()
    return samples, fastq_fetch_params, star_fastq_files, read_type

def check_if_ref_files_match(*, dna: FileName, rna: FileName, species: define_constants().known_species):

    """
    Input = paired URLs
    https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

    https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
    """

    species_Latin = {"human" : "sapiens", "mouse": "musculus"}.get(species)
    if species_Latin in dna and species_Latin in rna:
        log_message(f"DNA and RNA files for {species} both contain {species_Latin}")
    else:
        log_message(f"Mismatch in DNA and RNA files for {species} aka {species_Latin}:\nDNA: {dna}\nRNA{rna}", exit_now=True)

    dna_build = os.path.basename(dna).split(".")[1]
    dna_build_from_rna_file, rna_build = os.path.basename(rna).split(".")[1:3]

    if dna_build == dna_build_from_rna_file:
        log_message(f"DNA and RNA files are both for genome build {dna_build}")
    else:
        log_message(f"Genome builds mismatched {dna_build} vs. {dna_build_from_rna_file}:\nDNA: {dna}\nRNA{rna}", exit_now=True)
    
    ens = []
    for f in [dna, rna]:
        if m := re.search(r'release-(\d+)\/', f):
            ens.append(m.groups()[0])
        else:
            log_message(f"Release not found in {f}", exit_now=True)
    if len(set(ens)) != 1:
        log_message(f"Ensembl releases mismatched {ens[0]} vs. {ens[1]}:\nDNA: {dna}\nRNA{rna}", exit_now=True)
    ens = ens[0]
    return dna_build, ens
    

def star_make_idx(args):
    # Output command to make an index for a given species and read length.
    if args.SCRIPT:
        args.SCRIPT = args.SCRIPT.replace("{species}", args.SPECIES)
        args.SCRIPT = args.SCRIPT.replace("{readlength}", str(args.READLENGTH))
        if not args.OVERWRITE:
            exit_if_files_exist(args.SCRIPT)
    index = args.INDEX
    index = index.replace("{species}", args.SPECIES)
    index = index.replace("{readlength}", str(args.READLENGTH))

    dna = constants.ref.dna[args.SPECIES]
    rna = constants.ref.rna[args.SPECIES]
    print(f"dna={dna}, rna={rna}, species={species}")
    dna_build, rna_build = check_if_ref_files_match(dna=dna, rna=rna, species=args.SPECIES)
    #"human": "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz",
    #"mouse": "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
    dna = os.path.basename(dna).replace(".gz","")
    rna = os.path.basename(rna).replace(".gz","")
    index = index.replace("{dna_build}", dna_build)
    index = index.replace("{rna_build}", rna_build)

    cpus = args.CPUS or "$(grep -c ^processor /proc/cpuinfo)"
    if args.MEM == 0:
        mem=""
    else:
        mem = f"--limitGenomeGenerateRAM {args.MEM}"
    output = [bash_header()]
    output.append(dedent(f'''
        refdir={constants.ref_dir}
        dna=$refdir/{dna}
        rna=$refdir/{rna}
        index={index}
        cpus={cpus}
        if [[ -d $index ]]; then
            die "$index exists"
        fi
        for file in $dna $rna ; do
            if [[ ! -e $file ]]; then
                if [[ -e $file.gz ]]; then
                    gunzip --keep $file.gz
                else
                    die "$file not found"
                fi
            fi
        done
        STAR --runMode genomeGenerate --genomeDir $index \\
            --genomeFastaFiles $dna --sjdbGTFfile $rna --sjdbOverhang {args.READLENGTH} \\
            --runThreadN $cpus {mem} > $index.log 2> $index.err
    '''))
    if args.SCRIPT:
        with open(args.SCRIPT, "w") as tempout:
            print("\n".join(output), file=tempout)
        os.chmod(args.SCRIPT, 0o755)
        log_message(f"created {args.SCRIPT}")
    else:
        print("\n".join(output))

def star(args):
    # output commands to fetch fastqs and run STAR and featureCounts
    test_executables(["samtools", "STAR"], exit_on_error=True)
    verify_that_paths_exist(
        [constants.star.base_dir, constants.bin_dir, constants.ref_dir],
        exit_on_error=True,
    )
    cpus = args.CPUS
    sortcpus = args.SORT_CPUS
    sortmem = args.SORTMEM
    # prefetch = args.PREFETCH
    verify_that_paths_exist(args.INPUTFILE)
    if args.ABRA2COORDS:
        verify_that_paths_exist(args.ABRA2COORDS)
        test_executables("abra2", exit_on_error=True)
    if args.RNASPADESCOORDS:
        verify_that_paths_exist(args.RNASPADESCOORDS)
        test_executables(rnaspades.exe)

    outputdir = check_dir(dir=args.OUTPUTDIR)
    counts_dir = check_dir(dir=os.path.join(outputdir, "counts"))
    if glob.glob(f"{counts_dir}/*"):
        log_message(f"\nCaution: files found in {counts_dir}\n", exit_now=True)
    if args.BAMDESTDIR:
        bamdestdir = check_dir(dir=args.BAMDESTDIR)
    else:
        bamdestdir = ""
    if args.ABRA2COORDS:
        abra2dir = check_dir(dir=os.path.join(outputdir, "abra2"))
    if args.RNASPADESCOORDS:
        rnaspadesdir = check_dir(dir=os.path.join(outputdir, "rnaspades"))
        # rnaspadesdir = os.path.join(outputdir, "rnaspades")
        # temp = os.path.basename(os.path.realpath(outputdir))
        # rnaspadesdir = check_dir(dir = os.path.join(rnaspadesdir, temp))

    outputscript = args.SCRIPT  # os.path.join(outputdir, )
    overwrite = args.OVERWRITE
    if outputscript and not overwrite:
        exit_if_files_exist(outputscript)
    indexpath, species = find_star_index(args.INDEX)
    # find gtf
    gtf = find_gtf(species=species)
    fastq_source = args.FASTQSOURCE.upper()
    samples, fastq_fetch_params, star_fastq_files, read_type = parse_fastq_manifest(
        inputfile=args.INPUTFILES,
        fastq_source=fastq_source,
        outputdir=outputdir,
        overwrite=overwrite,
    )

    tempbam = {
        sample: os.path.join(outputdir, f"{sample}.Aligned.out.bam")
        for sample in samples
    }
    exit_if_files_exist(list(tempbam.values()))
    sortedbam = {sample: os.path.join(outputdir, f"{sample}.bam") for sample in samples}
    exit_if_files_exist(list(sortedbam.values()))

    if bamdestdir:
        bamfiles = [
            os.path.join(bamdestdir, os.path.basename(x)) for x in sortedbam.values()
        ]
        bamindexes = [f"{x}.bai" for x in bamfiles]
        exit_if_files_exist(bamfiles + bamindexes)
    countsfile = {
        sample: os.path.join(counts_dir, f"{sample}.counts.txt") for sample in samples
    }
    exit_if_files_exist(list(countsfile.values()))

    output = []
    output.append(bash_header())
    bash_functions = {}
    output.append(dedent(f"""
        # option to load index while fetching first set of fastqs
        # {constants.bin_dir}/load.index.into.mem.sh {indexpath} &>> load.index.log &
    """))

    if fastq_source == "ENA":
        function_name, function_code = bash_aspera(speed=args.TRANSFER_SPEED)
        bash_functions["getfastqs"] = function_name
        output.append(function_code)
        readFilesCommand = "zcat"
    else:  # SRA
        function_name, function_code = bash_fasterq(
            cpus=int(cpus / 2), read_type=read_type
        )
        bash_functions["getfastqs"] = function_name
        output.append(function_code)
        readFilesCommand = "cat"

    function_name, function_code = bash_star(
        cpus=args.CPUS, index=indexpath, readFilesCommand=readFilesCommand
    )
    bash_functions["star"] = function_name
    output.append(function_code)

    function_name, function_code = bash_samtools_sort_by_position(
        sortcpus=sortcpus, sortmem=sortmem
    )
    bash_functions["sort BAM by read ID"] = function_name
    output.append(function_code)

    if args.COUNTS:
        function_name, function_code = bash_featurecounts(
            gtf=gtf, cpus=cpus, options=constants.featurecounts_options
        )
        bash_functions["featurecounts"] = function_name
        output.append(function_code)

    for sample in samples:
        output.append(f"\n#{sample}")
        if bamdestdir:
            temp = os.path.join(bamdestdir, os.path.basename(sortedbam[sample]))
            output.append(dedent(f"""
                # bamdest set to {bamdestdir} - check for output bam {temp}
                if [[ -e {temp} ]] ; then 
                    die "conflict - found {temp}"
                fi
            """))
        getfastqs_cmd = bash_functions["getfastqs"]
        getfastqs_params = fastq_fetch_params[sample]
        outputprefix = os.path.join(outputdir, f"{sample}")
        star_cmd = bash_functions["star"]
        star_params = f'"{star_fastq_files[sample]}" {outputprefix}'
        other_opts = []
        if args.ADD_OPTS:
            other_opts.append(args.ADD_OPTS)
        if args.ADDREADGROUP:
            other_opts.append(f'--outSAMattrRGline \\"ID:{sample}\\"')
        other_opts = " ".join(other_opts)
        star_params += f' "{other_opts}"'
        files_to_delete = star_fastq_files[sample].replace(",", " ")
        output.append(dedent(f"""
            # get fastqs
            {getfastqs_cmd} {getfastqs_params}
            # 
            # run  STAR
            {star_cmd} {star_params}
            #
            # cleanup
            rm {files_to_delete}
            #
        """))
        if args.COUNTS:
            featurecounts_cmd = bash_functions["featurecounts"]
            featurecounts_params = " ".join([tempbam[sample], countsfile[sample]])
            output.append(dedent(f"""
                # counts
                {featurecounts_cmd} {featurecounts_params}
            """))

        output.append("\n# sort BAM and index")
        sort_bam_cmd = bash_functions["sort BAM by read ID"]
        sort_bam_params = " ".join([tempbam[sample], sortedbam[sample]])
        output.append(dedent(f"""
            # sort the BAM file by position
            {sort_bam_cmd} {sort_bam_params}
            samtools index -@ {cpus} {sortedbam[sample]} &>> {outputprefix}.index.log
            #
            # delete the original
            rm {tempbam[sample]}
        """))
        if bamdestdir:
            output.append(dedent(f"""
                # move files in background task
                nice mv --no-clobber {sortedbam[sample]} {sortedbam[sample]}.bai {bamdestdir}/ &>> {outputprefix}.mv.log &
                # keep track of the job  ID
                mvjobid=$!
        """))
        output.append(dedent(f"""
            zip -qymT star.outputs.and.logs.zip {sample}.SJ.out.tab {sample}.Log.final.out {sample}.Log.out {sample}.Log.progress.out {sample}.star.out &>> {outputprefix}.zip.log
        """))

    find_err_files = f"ls {outputdir}*.err {outputdir}*.log 2> /dev/null"
   
    output.append(dedent(f"""
        set +e
        # clear mem
        {constants.bin_dir}/remove.genome.from.mem.sh {indexpath} &
        sed -ri '/^#|bam_sort_core..merging|nohup..ignoring.input/d' $({find_err_files}) &>> /dev/null
        find {outputdir} -type f -empty -delete
        {__file__} counts_postproc --species {species} --dir {counts_dir} --exec
        {__file__} star_zip_files
        {__file__} star_clear_mem --index {indexpath}
    """))

    if bamdestdir and (args.ABRA2COORDS or args.RNASPADESCOORDS):
        output.append(dedent(f"""
            log waiting for mvjobid=$mvjobid
            wait $mvjobid
            log $mvjobid done
            # temp = os.path.realpath(bamdestdir)
        """))

    if args.ABRA2COORDS:
        srcdir = bamdestdir or os.path.realpath(outputdir)
        output.append(f"cp -s {srcdir}/*.bam {srcdir}/*.bai {abra2dir}/")
        output.append(f"mv {args.ABRA2COORDS} {abra2dir}/")
        output.append(f"cd {abra2dir}")
        output.append(
            f"{__file__} abra2 -b *.bam -r {os.path.basename(args.ABRA2COORDS)} -z --exec"
        )
        output.append(f"cd {os.path.realpath(outputdir)}")

    if args.RNASPADESCOORDS:
        srcdir = bamdestdir or os.path.realpath(outputdir)
        output.append(f"cp -s {srcdir}/*.bam {srcdir}/*.bai {rnaspadesdir}")
        output.append(f"mv {args.RNASPADESCOORDS} {rnaspadesdir}")
        output.append(f"cd {rnaspadesdir}")
        output.append(
            f"{__file__} rnaspades -b *.bam -r {os.path.basename(args.RNASPADESCOORDS)} -z --exec"
        )
        output.append(f"cd {os.path.realpath(outputdir)}")

    cmd = f"nohup {outputscript} &>> {outputscript}.log &" if outputscript else ""
    handle_commands_and_output(
        commands=output, outputfile=outputscript, single_line_command=cmd, execute =False
    )


def featurecounts_ids(args):
    """
    Determines mapping to unique column names when inputs are combined
    """
    files = args.INPUTFILES
    verify_that_paths_exist(files)
    outputfile = args.OUTPUTFILE
    if outputfile:
        outputdir = args.OUTPUTDIR.rstrip("/")
        outputfile = os.path.join(outputdir, outputfile)
        exit_if_files_exist(outputfile)
    output = []
    newsamples = []
    output.append("sample\tnew_sample\tfile")
    for file in files:
        data = pd.read_csv(file, sep="\t", header=0, comment="#", nrows=1)
        for sample in set(data.columns) - set(constants.featurecounts_constant_columns):
            newsample = sample
            for ext in [".Aligned.out.bam", ".unsorted.bam", ".bam"]:
                if newsample.endswith(ext):
                    newsample = re.sub(ext, "", newsample)
            newsample = newsample.split("/")[-1]
            output.append(f"{sample}\t{newsample}\t{file}")
            newsamples.append(newsample)

    if len(newsamples) != len(set(newsamples)):
        log_message("\nCaution: non-unique sample IDs\n", exit_now=True)

    if outputfile:
        with open(outputfile, "w") as tempout:
            print("\n".join(output), file=tempout)
        log_message(f"Created {outputfile}")
    else:
        print("\n".join(output))


def get_one_file_of_counts(file=None):
    # gets geneid and counts from featurecounts - skips all metadata
    data = pd.read_csv(file, sep="\t", header=0, comment="#", nrows=0)
    columns_to_get = list(
        set(data.columns) - set(constants.featurecounts_constant_columns[1:])
    )
    data = pd.read_csv(file, sep="\t", header=0, comment="#", usecols=columns_to_get)
    return data

def sum_counts(data=None, sum_column="sum", outputfile=None, overwrite=False):

    # input = dataframe
    # returns sum by sample as an OrderedDict
    # optionally outputs sums to a file
    # if  outputfile is None, assumes stdout
    if outputfile and overwrite is False and os.path.exists(outputfile):
        log_message(
            f"{outputfile} exists - delete it or allow overwrite", exit_now=True
        )
    # assert outputfile is not None, "malformed outputfile"
    samples = [
        x for x in data.columns[1:] if not x in constants.featurecounts_constant_columns
    ]
    total_counts = OrderedDict({sample: np.sum(data[sample]) for sample in samples})
    with open(outputfile, "w") if outputfile else sys.stdout as tempout:
        print(f"sample\t{sum_column}", file=tempout)
        for sample in samples:
            print(f"{sample}\t{int(total_counts[sample])}", file=tempout)
    return total_counts


def sumcounts(args):
    log_message("sum counts")
    verify_that_paths_exist(args.INPUTFILE)
    outputfile = args.OUTPUTFILE
    # exit_if_files_exist(outputfile)
    # header = None if args.NOHEADER else 0
    data = pd.read_csv(args.INPUTFILE, sep="\t", header=0)
    sum_counts(
        data=data,
        sum_column="total_counts",
        outputfile=args.OUTPUTFILE,
        overwrite=args.OVERWRITE,
    )
    log_message(f"Total counts written to {outputfile or 'stdout'}")


def add_gene_name(args):
    log_message("substitute gene name in expression counts if unique")
    verify_that_paths_exist(args.INPUTFILE)
    outputfile = args.OUTPUTFILE
    if outputfile:
        exit_if_files_exist(outputfile)
    species = args.SPECIES
    ref_file = f"${HOME}/star/ref/ens.to.gene.{species}.if.unique"
    gene_names = pd.read_csv(
        ref_file, sep="\t", header=None, names=["Geneid", "dummy"], dtype=object
    )
    data = pd.read_csv(args.INPUTFILE, sep="\t", header=0, dtype=object)
    data = pd.merge(data, gene_names, on="Geneid", how="inner")
    data.insert(0, "gene_name", "")
    data["gene_name"] = data["dummy"]
    data.drop(columns=["Geneid", "dummy"]).to_csv(
        outputfile or sys.stdout, sep="\t", index=False
    )
    temp = outputfile or "stdout"
    log_message(f"counts by gene name written to {temp}")


def merge_counts(args):
    log_message("merge counts")
    files = args.INPUTFILES
    verify_that_paths_exist(files)
    outputfile = args.OUTPUTFILE
    outputdir = args.OUTPUTDIR.rstrip("/")
    totals = args.TOTALS
    if outputdir:
        if outputfile and os.sep in outputfile:
            log_message(
                "outputfile cannot include dir if outputdir is specified", exit_now=True
            )
        if totals and os.sep in totals:
            log_message(
                "output file for totals cannot include dir if outputdir is specified",
                exit_now=True,
            )
        totals = args.TOTALS
        outputfile = os.path.join(outputdir, outputfile)
        totals = os.path.join(outputdir, totals)
        exit_if_files_exist([outputfile, totals])
    output = get_one_file_of_counts(file=files[0])
    for file in files[1:]:
        data = get_one_file_of_counts(file=file)
        if common_samples := set(data.columns[1:]) & set(output.columns[1:]):
            common_samples = "\n".join(common_samples)
            log_message(f"Common sample IDs:\n{common_samples}", exit_now=True)
        output = output.merge(data, how="inner", on="Geneid")
    rename_samples = args.RENAME_SAMPLES
    if rename_samples:
        rename_samples = pd.read_csv(
            args.RENAME_SAMPLES,
            header=None,
            sep="\t",
            usecols=[0, 1],
            index_col=0,
        ).T.to_dict("records")[0]
        output.rename(columns=rename_samples, inplace=True)
    output.to_csv(outputfile, sep="\t", index=False)
    temp = outputfile or "stdout"
    log_message(f"Combined counts written to {temp}")
    sum_counts(data=output, sum_column="total_counts", outputfile=totals)
    log_message(f"Total counts written to {totals}")


def exp(args):
    log_message("transform counts")
    inputfile = args.INPUTFILE


def transformcounts(args):
    log_message("transform counts")
    verify_that_paths_exist(args.INPUTFILE)
    method = args.METHOD.lower()
    if "log2" in method and args.LOG2_OFFSET <= 0:
        log_message(
            "log2_offset must be positive when applying log", exit_now=True
        )

    outputfile = args.OUTPUTFILE
    if not outputfile:
        outputfile = f"{method.lower()}.txt"
    overwrite = args.OVERWRITE
    if outputfile and overwrite is False:
        exit_if_files_exist(outputfile)

    # get metadata
    metadata = args.METADATA
    if metadata:
        metadata = os.path.join(os.path.dirname(args.INPUTFILE), metadata)
        verify_that_paths_exist(metadata)
        use_genenames = args.GENE_NAMES if metadata else 0
        metadata_columns_to_get = ["Geneid", "Length", "gene_biotype"]
        if use_genenames:
            metadata_columns_to_get.append("gene_name")
        temp = pd.read_csv(metadata, sep="\t", header=0, comment="#", nrows=1)
        verify_columns_present_in_dataframe(
            data=temp, columns=metadata_columns_to_get, source=metadata
        )
        metadata = pd.read_csv(
            metadata, sep="\t", header=0, usecols=metadata_columns_to_get, comment="#"
        )
        gene_types = args.GENE_TYPES
        if gene_types and gene_types != "all":
            mask = metadata["gene_biotype"].isin(gene_types)
            metadata = metadata[mask].copy()
    data = pd.read_csv(args.INPUTFILE, sep="\t", header=0)
    #  sum counts in case gene IDs are repeaed
    data = data.groupby(data.columns[0], as_index=False).sum()
    samples = [
        x for x in data.columns[1:] if not x in constants.featurecounts_constant_columns
    ]

    # ["CPM", "CPM-UQ", "CPM-UQ-log2", "RPKM", "RPKM-UQ", "RPKM-UQ-log2"]
    temp_output = f"{args.TOTALS}.txt"
    total_counts = sum_counts(
        data=data,
        sum_column="total_counts",
        outputfile=temp_output,
        overwrite=overwrite,
    )

    if metadata:
        before = data.shape[0]
        data = data.merge(metadata, how="inner", on="Geneid")
        log_message(
            f"{data.shape[0]} of {before} genes remain after merging data & metadata"
        )
        if use_genenames:
            mask = data["gene_name"].isnull()
            data.loc[~mask, "Geneid"] = data["gene_name"]
            data.drop(columns=["gene_biotype", "gene_name"], inplace=True)
        else:
            data.drop(columns=["gene_biotype"], inplace=True)
        if use_genenames:
            data = data.groupby(data.columns[0], as_index=False).sum()

        temp_output = f"{args.TOTALS}_filtered.txt"
        total_counts = sum_counts(
            data=data, sum_column="filtered_totals", outputfile=temp_output
        )

    samples = [
        x for x in data.columns[1:] if not x in constants.featurecounts_constant_columns
    ]
    if method == "percentile":
        for col in samples:
            data[col] = data[col].rank(pct=True) * 100
            data[col] = data[col].round(0).astype(int)
    else:
        for sample in samples:
            if method.startswith("rpkm"):
                data[sample] *= 1000000000 / data["Length"] / total_counts[sample]
            elif method.startswith("cpm"):
                data[sample] *= 1000000 / total_counts[sample]
            else:
                log_message(f"Unknown method {method}", exit_now=True)

        temp_output = f"{args.TOTALS}_sum_{method}.txt"
        sum_counts(
            data=data,
            sum_column=f"sum_{method}",
            outputfile=temp_output,
            overwrite=overwrite,
        )

    if metadata:
        if columns_to_drop := set(["Length", "gene_biotype"]) & set(data.columns):
            data.drop(columns=list(columns_to_drop), inplace=True)
    if "uq" in method in method:
        log_message("Calculating upper quantile.")
        tempdata = (
            data[["Geneid"] + samples]
            .set_index("Geneid")
            .applymap(lambda x: x if x > 0 else np.NAN)
        )
        pct = np.nanpercentile(tempdata, args.RESCALE_PERCENTILE, axis=0)
        data[samples] = data[samples] * args.RESCALE_COMMON_VALUE / pct

    if "log2" in method:
        data[samples] = np.log2(data[samples].applymap(lambda x: x + args.LOG2_OFFSET))

    data.round(2).to_csv(outputfile or sys.stdout, sep="\t", index=False)

    temp = outputfile or "stdout"
    log_message(f"{method} written to {temp}")


def counts_postproc(args):
    output = []
    if args.SCRIPT:
        exit_if_files_exist(args.SCRIPT)
    dir = args.DIR.rstrip("/")
    species = args.SPECIES
    output = []
    output.append(bash_header())
    output.append(dedent(f"""
        script={__file__}
        if [[ -e {dir}/counts.txt ]] ; then
            input="{dir}/counts.txt"
            metadata_input=$input
            zip="zip -qymT {dir}/counts.zip {dir}/counts.txt {dir}/counts.txt.summary {dir}/counts.log"
        elif [[ $(ls {dir}/*.counts.txt 2> /dev/null) ]] ; then
            input="{dir}/*.counts.txt"
            metadata_input=$(ls {dir}/*.counts.txt | head -n 1)
            zip="zip -qymT {dir}/counts.zip {dir}/*.counts.*"
        else
            die "No counts files found in {dir}"
        fi
        
        # resolve unique mapping of samples before merging files
        $script featurecounts_ids -i $input --outputdir {dir} 
        
        # merge counts results for samples from this experiment
        $script merge_counts -i $input --rename_samples {dir}/{constants.counts_sample_ID_file} --outputdir {dir} --outputfile {constants.default_merged_counts}
        
        # add the gene name to Ensembl IDs if one-to-one
        $script add_gene_name --inputfile {dir}/{constants.default_merged_counts} --outputfile {dir}/{constants.default_counts_by_gene_name} --species {species}

        # extract gene metadata from one of the files 
        $script gene_metadata -i $metadata_input --outputdir {dir} --outputfile {constants.default_gene_metadata}
        
        # zip individual files after merging
        $zip
        gzip {dir}/{constants.default_merged_counts} {dir}/{constants.default_gene_metadata}
    """))

    # output = ["#!/bin/bash\n"] + output
    cmd = f"{args.SCRIPT} &>> {args.SCRIPT}.log" if args.SCRIPT else ""
    handle_commands_and_output(
        commands=output,
        outputfile=args.SCRIPT,
        single_line_command=cmd,
        execute =args.EXEC,
    )


def star_zip_files(args):
    dir = args.DIR.rstrip("/")
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    # find_err_files = f"find {dir} -maxdepth 1 -name \"*.err\""

    if dir == ".":
        dir = ""
    else:
        dir += "/"
    output = []
    output.append(bash_header(flags="")) #"#!/bin/bash")  # avoid set -e

    # find_err_files = f"find {dir} -maxdepth 1 -name \"*.err\" "
    find_err_files = f"ls {dir}*.err {dir}*.log 2> /dev/null"
    output.append(dedent(f"""
        sed -ri '/bam_sort_core..merging|nohup..ignoring.input/d' $({find_err_files} 2> /dev/null
        {find_err_files} -empty -delete
        zip -qymT {dir}star.outputs.and.logs.zip {dir}*SJ*tab {dir}*Log* {dir}*.out
        ls {dir}*star*sh {dir}*manifest* {dir}*fastq_file* 2> /dev/null |grep -v zip$ |zip -qymT -@ star.commands.and.manifests.zip
    """))
    # output.append(f"ls fastq.* 2> /dev/null |grep -v zip$ |zip -T -@ star.commands.and.manifests.zip")
    cmd = (
        f"{outputscript} &>> {outputscript}.log 2> {outputscript}.err"
        if outputscript
        else ""
    )
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )

def star_clear_mem(args):
    outputscript = args.SCRIPT
    if outputscript:
        exit_if_files_exist(outputscript)
    verify_that_paths_exist(constants.star.dummy_fastq)
    output = [bash_header()]
    output.append(dedent(f"""
        STAR --genomeDir {args.INDEX} --readFilesIn {constants.star.dummy_fastq} --runThreadN 4 --outFileNamePrefix ${constants.star.dummy_fastq}. --genomeLoad Remove --runMode alignReads"
    """))
    # --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 25 --seedMultimapNmax 250 --seedPerReadNmax 250 --alignTranscriptsPerReadNmax 1000 --limitOutSJoneRead 100 --alignWindowsPerReadNmax 250 
    output = output.dedent()
    with open(outputscript, "w") if outputscript else sys.stdout as tempout:
        print(output, file=tempout)
    if outputscript:
        os.chmod(outputscript, 0o755)
        log_message(f"created {outputscript}")
        #print(outputscript)

def check_for_file_xor_stdin(inputfile: FileName | None):
    if inputfile is None:
        if detect_stdin():
            return
        log_message("Specify an input file or stdin.", exit_now=True)
    elif detect_stdin():
        log_message("Specify an input file or stdin, but not both.", exit_now=True)
        verify_that_paths_exist(args.INPUTFILE)


def text_to_fasta(args):
    # inputs tab-delimited text, outputs fasta
    check_for_file_xor_stdin(args.INPUTFILE)
    outputfile = args.OUTPUTFILE
    if outputfile:
        exit_if_files_exist(outputfile)
    # seq = ""
    minlength = args.MIN_LENGTH
    with open(args.INPUTFILE, "r") if args.INPUTFILE else sys.stdin as tempin:
        with open(outputfile, "w") if outputfile else sys.stdout as tempout:
            line = tempin.readline()
            temp = line.rstrip().split("\t")
            if temp[1] != "sequence":
                if minlength == 0 or len(temp[1]) >= minlength:
                    print(f">{temp[0]}\n{temp[1]}", file=tempout)
            if minlength:
                for line in tempin:
                    temp = line.rstrip().split("\t")
                    if len(temp[1]) >= minlength:
                        # tempout.write(f">{temp[0]}\n{temp[1]}\n")
                        print(f">{temp[0]}\n{temp[1]}", file=tempout)
            else:
                for line in tempin:
                    line = re.sub("\t", "\n", line)
                    print(f">{line.rstrip()}", file=tempout)
                    # tempout.write(">" + line)


def fasta_to_text(args):
    # inputs fasta file, outputs text file with header row
    # presumably called from command line
    # output to a file
    if args.INPUTFILE:
        verify_that_paths_exist(args.INPUTFILE)
    else:
        detect_stdin(exit_on_error=True)
    outputfile = args.OUTPUTFILE
    if outputfile:
        exit_if_files_exist(outputfile)
    minlength = args.MIN_LENGTH
    seq = ""
    if minlength:
        log_message(f"Caution: min length = {minlength}")
    with open(args.INPUTFILE, "r") if args.INPUTFILE else sys.stdin as tempin:
        with open(outputfile, "w") if outputfile else sys.stdout as tempout:
            print("ID\tsequence", file=tempout)
            for line in tempin:
                line = line.rstrip()
                if line.startswith(">"):
                    if seq and len(seq) > minlength:
                        print(f"{id}\t{seq}", file=tempout)
                        # tempout.write("%s\t%s\n" % (id,seq))
                    id = line[1:]
                    seq = ""
                else:
                    seq += line
            if seq and len(seq) > minlength:
                # tempout.write("%s\t%s\n" % (id,seq))
                print(f"{id}\t{seq}", file=tempout)


'''
def bash_check_transfer_speed(file_name="aspera.speed"):
    # use "aspera.speed" to allow changes in transfer speed
    bash_function_name = "update_transfer_speed"
    code = f"""
        {bash_function_name}(){
            file={file_name}
            if [[ -e $file ]] ; then speed=$(cat $file)m ; echo $speed ; fi
        }"""
    code = dedent(text=code)
    return bash_function_name, code
'''


def dedup_cols(data):
    assert isinstance(data, pd.DataFrame), "data is not a dataframe"
    column_contents = {x: "\n".join(map(str, data[x])) for x in data.columns}
    skip = {x: 0 for x in data.columns}
    rename = {}
    for i, col1 in enumerate(data.columns[0:-1]):
        if skip[col1]:
            continue
        for col2 in data.columns[i + 1 :]:
            if skip[col2]:
                continue
            if column_contents[col1] == column_contents[col2]:
                skip[col2] = 1
                # rename columns that came from merging if dropping one
                if (
                    col2.endswith("_y")
                    and col1.endswith("_x")
                    and col1[0:-2] == col2[0:-2]
                ):
                    rename[col1] = col1[0:-2]
                log_message(f"{col2} repeats {col1}")
    if cols_to_drop := [x for x in data.columns if skip[x]]:
        data.drop(columns=cols_to_drop, inplace=True)
    if cols_to_rename := set(rename.keys()) & set(data.columns):
        rename = {x: rename[x] for x in cols_to_rename}
        if len(rename.keys()) != len(rename.values()):
            log_message("Column name conflict")
        else:
            data.rename(columns=rename, inplace=True)
    return data


def concat_files(args):
    if len(args.INPUTFILES) < 2:
        log_message("Specify multiple input files", exit_now=True)
    verify_that_paths_exist(args.INPUTFILES)
    outputfile = args.OUTPUTFILE
    if outputfile:
        exit_if_files_exist(outputfile)
    data = pd.read_csv(args.INPUTFILES[0], sep="\t", header=0, dtype=object)
    for file in args.INPUTFILES[1:]:
        temp = pd.read_csv(file, sep="\t", header=0, dtype=object)
        data = pd.concat([data, temp], sort=False, ignore_index=True)
    data.to_csv(outputfile if outputfile else sys.stdout, sep="\t", index=False)


def join_files(args):
    Nfiles = len(args.INPUTFILES)
    if Nfiles < 2:
        log_message("Specify multiple input files", exit_now=True)
    verify_that_paths_exist(args.INPUTFILES)
    outputfile = args.OUTPUTFILE
    if outputfile:
        exit_if_files_exist(outputfile)
    method = args.METHOD
    columns = args.COLUMNS
    if len(columns) == 1:
        # use same columns for all files
        cols = columns[0].split(",")
        columns = [cols for i in range(Nfiles)]
    elif len(columns) == Nfiles:
        for i, cols in enumerate(columns):
            columns[i] = cols.split(",")
    else:
        log_message(
            f"Invalid # columns for {Nfiles} files. Specify either one or {Nfiles} column or list of columns",
            exit_now=True,
        )
    for i, inputfile in enumerate(args.INPUTFILES):
        data = pd.read_csv(inputfile, sep="\t", header=0, dtype=object)
        # data[i].to_csv(f"out.{i}", sep="\t", index=False)
        # print(f"{i}:\n" + "\n".join(data[i].columns))
        verify_columns_present_in_dataframe(
            data=data, columns=columns[i], source=inputfile
        )
        if i:
            output = pd.merge(
                output, data, left_on=columns[0], right_on=columns[i], how=method
            )
        else:
            output = data
        if output.empty:
            log_message(
                f"No common column values after reading {inputfile}", exit_now=True
            )
    if args.DEDUP:
        output = dedup_cols(data=output)
    output.to_csv(outputfile if outputfile else sys.stdout, sep="\t", index=False)


"""
def join(args):
    inputfiles = args.INPUTFILES
    if len(inputfiles) != 2:
        log_message("Specify two input files", exit_now=True)
    verify_that_paths_exist(inputfiles)
    outputfile = args.OUTPUTFILE
    if outputfile: exit_if_files_exist(outputfile)
    method = args.METHOD
    columns = args.COLUMNS
    if len(columns) == 2:
        for i, cols in enumerate(columns):
            columns[i] = cols.split(",")
    elif len(columns) == 1:
        columns.append(columns[0].split(","))
        columns[0] = columns[0].split(",")
        #columns += columns[0]
    else:
        log_message("Invalid columns", exit_now=True)
    data = []
    for i in [0, 1]:
        data.append(pd.read_csv(inputfiles[i], sep="\t", header = 0, dtype=object))
        #data[i].to_csv(f"out.{i}", sep="\t", index=False)
        #print(f"{i}:\n" + "\n".join(data[i].columns))
        verify_columns_present_in_dataframe(data = data[i], columns=columns[i], source = inputfiles[i])
    output = pd.merge(data[0], data[1], left_on = columns[0], right_on=columns[1], how = method)
    if output.empty:
        log_message("No common column values", exit_now=True)
    if args.DEDUP:
        output = dedup_cols(data=output)
    output.to_csv(outputfile if outputfile else sys.stdout, sep="\t", index=False)
"""


def minimap_output_file(inputfile: FileName):  # =None, ext=None):
    outputfile = re.sub(r"\.gz$", "", inputfile)
    outputfile = re.sub(".fasta", "", outputfile)
    # outputfile = re.sub(r"\.fa$", "", outputfile)
    # if ext:
    #    outputfile = f'{outputfile}.{ext.lstrip(".")}'
    return outputfile


def bash_minimap(*, cpus: int, bases_loaded: str, genome: str, junctions: FileName):
    # assert bases_loaded is not None, "Specify bases_loaded"
    # assert cpus is not None, "Specify cpus"
    # -2 Use two I/O threads during mapping.
    # -a : outputformat = sam
    bash_function_name = "run_minimap"
    command = "minimap2"
    # [ERROR] dual gap penalties violating E1>E2 and O1+E1<O2+E2
    # Gap extension penalty [2,1]. A gap of length k costs min{O1+k*E1,O2+k*E2}. In the splice mode, the second gap penalties are not used.
    # options = "-a -c --eqx --MD --MD --sr -xsplice -ub -G500k --junc-bonus=5 -C3 -O2,4 -E1,0" # -p 0.95 --secondary=no
    # options = "--secondary=no --eqx --MD -ax splice -G500k --junc-bonus=5 -C3"  # --sr -O2,4 -E1,0
    options = "--secondary=no --eqx --MD -ax splice -G500k --junc-bonus=5 -C3 -O2,4 -E1,0"  # --sam-hit-only"

    if bases_loaded:
        assert bases_loaded[-1].lower() in [
            "k",
            "m",
            "g",
        ], f"bases_loaded ends with {bases_loaded[-1]}, expected k, m or g"
        assert bases_loaded[
            :-1
        ].isdigit(), f"invalid bases_loaded - expecting an integer instead of {bases_loaded[:-1]}"
        options += f" -K{bases_loaded}"
    if cpus:
        cpus_sam = f"-t {cpus}"
        cpus_bam = f"-t {cpus - 1}"
        # assert str(cpus).isdigit(), f"invalid cpus {cpus}"
        # cpus -= 1 # account for samtools view  & sort
        # command += f" -t{cpus}"
    else:
        cpus_sam = cpus_bam = ""
    code = f"""
        {bash_function_name}(){{
            inputfile=$1
            outputfileprefix=$2
            outputformat=$3
            genome={genome}
            junctions="--junc-bed {junctions}"
            options="{options}"
            if [[ $outputformat == "sam" ]] ; then
                outputfile=$outputfileprefix.sam
                if [[ -e $outputfile ]] ; then log "$outputfile exists" ; return ; fi
                minimap2 $options {cpus_sam} $junctions -o $outputfile $genome $inputfile 2> $outputfile.log
            elif [[ $outputformat == "bam" ]] ; then
                outputfile=$outputfileprefix.bam
                if [[ -e $outputfile ]] ; then log "$outputfile exists" ; return ; fi
                minimap2 $options {cpus_bam} $junctions $genome $inputfile 2> $outputfile.1.log | samtools view -h -q 1 | samtools sort -o $outputfile - 2> $outputfile.2.log
                samtools index $outputfile
                if type "bam2Bed12" ; then
                    bam2Bed12 -i $outputfile > $outputfile.bed
                fi
            else
                die "unknown output format $outputformat"
            fi
        }}"""
    code = dedent(text=code)
    return bash_function_name, code


def minimap(args):
    # needs coords and either gene bams or non-gene BAMs, but not both
    test_executables("minimap2")
    # outputfileext = args.EXT
    verify_that_paths_exist(args.INPUTFILES)
    outputscript = args.SCRIPT
    exit_if_files_exist(outputscript)
    species = args.SPECIES
    junctions = constants.minimap_junctions[species]
    genome = constants.minimap_genome[species]
    if genome == "":
        log_message(f"No minimap genome index for {species}", exit_now=True)
    if junctions == "":
        log_message(f"No minimap junctions for {species}", exit_now=True)
    genome = os.path.join(constants.ref_dir, genome)
    junctions = os.path.join(constants.ref_dir, junctions)
    verify_that_paths_exist([genome, junctions])
    outputformat = args.OUTPUTFORMAT
    # outputfileprefix = {x: minimap_output_file(inputfile=x, ext=outputfileext) for x in inputfiles}
    outputfileprefix = {x: minimap_output_file(inputfile=x) for x in args.INPUTFILES}
    exit_if_files_exist([f"{x}.{outputformat}" for x in outputfileprefix.values()])
    output = []
    output.append(bash_header())
    function_name, function_code = bash_minimap(
        cpus=args.CPUS, bases_loaded=args.MEM, genome=genome, junctions=junctions
    )
    output.append(function_code)
    for inputfile in args.INPUTFILES:
        output.append(f"# {inputfile}")
        output.append(
            f"{function_name} {inputfile} {outputfileprefix[inputfile]} {outputformat}"
        )
    cmd = f"nohup {outputscript} &>> {outputscript}.log &"
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.EXEC,
    )

    # bam2Bed12.py -i transcripts.fasta.v4.bam > transcripts.fasta.v4.bed2.txt

def ens_to_gene(args):
    gtf_chr =  0
    gtf_feature_type =  2
    gtf_features =  -1
    if args.SPECIES:
        if args.INPUTFILE or args.OUTPUTFILE:
            log_message("Specify species, or specify input and output files.", exit_now=True)
        inputfile = os.path.expandvars(constants.gtf[args.SPECIES])
        if not os.path.exists(inputfile):
            if os.path.exists(inputfile + ".gz"):
                cmd = f"gunzip --keep {inputfile}.gz"
                log_message(f"Executing {cmd}")
                execute_command(cmd)
        outputfile = constants.featurecounts.unique_IDs[args.SPECIES]
        outputfile = os.path.expandvars(outputfile)
    else:
        if args.INPUTFILE is None or args.OUTPUTFILE is None:
            log_message("Specify both input file and output file, or specify species.", exit_now=True)
        inputfile = args.INPUTFILE
        outputfile = args.OUTPUTFILE
    verify_that_paths_exist(inputfile)
    if not args.OVERWRITE:
        exit_if_files_exist(outputfile)
    # (inputfile: FileName, *, outputfile: FileName|None=None, sort: bool = False, subset=["protein_coding"]):
    gtf_features =  -1
    data = []
    with open(inputfile, "r") as tempin:
        for line in tempin:
            if line.startswith("#"):
                continue # skip comments
            #if not re.match(r'\d|chr\d|[MXY]|chr[MXY]', line):
            #    continue
            columns = line.rstrip().split("\t")
            if columns[gtf_feature_type] != "gene":
                continue
            features = columns[-1]
            output = []
            skip = False
            for kw in ['gene "|gene_name "', 'gene_id "']: #, 'gene_biotype "']:
                temp = re.split(kw, features)
                if len(temp) > 1:
                    output.append(temp[1].split('"')[0])
                else:
                    skip = True
            if skip:
                continue
            data.append(output)
    data = pd.DataFrame(data, columns = "gene_name gene_id".split()).sort_values(by="gene_name") #  gene_biotype
    ctr = Counter(data["gene_name"])
    unique = {k for (k, v) in ctr.items() if v == 1}
    mask  = data["gene_name"].isin(unique)
    data[mask][["gene_id", "gene_name"]].to_csv(outputfile, sep="\t", index=False)
    log_message(f"created {outputfile}")

def gtf_to_coords(args):
    if args.SPECIES:
        if args.INPUTFILE or args.OUTPUTFILE:
            log_message("Specify species, or specify input and output files.", exit_now=True)
        inputfile = os.path.expandvars(constants.gtf[args.SPECIES])
        if not os.path.exists(inputfile):
            if os.path.exists(inputfile + ".gz"):
                cmd = f"gunzip --keep {inputfile}.gz"
                log_message(f"Executing {cmd}")
                execute_command(cmd)
        outputfile = constants.coords_source_file[args.SPECIES]
        outputfile = os.path.expandvars(outputfile)
    else:
        if args.INPUTFILE is None or args.OUTPUTFILE is None:
            log_message("Specify both input file and output file, or specify species.", exit_now=True)
        inputfile = args.INPUTFILE
        outputfile = args.OUTPUTFILE
    verify_that_paths_exist(inputfile)
    if not args.OVERWRITE:
        exit_if_files_exist(outputfile)

    # (inputfile: FileName, *, outputfile: FileName|None=None, sort: bool = False, subset=["protein_coding"]):
    gtf_chr =  0
    gtf_start =  3
    gtf_stop =  4
    gtf_feature_type =  2
    data = []
    genes_skipped = set()
    chrs_skipped = set()
    with open(inputfile, "r") as tempin:
        for line in tempin:
            if line.startswith("#"):
                continue # skip comments
            #if re.match(r'#|KI|GL|JH', line): 
            if not re.match(r'\d|chr\d|[MXY]|chr[MXY]', line):
                chrs_skipped.add(line.split("\t", 1)[0])
                continue
            columns = line.rstrip().split("\t")
            if columns[gtf_feature_type] != "gene":
                continue
            features = columns[-1]
            output = [columns[i] for i in [gtf_chr, gtf_start, gtf_stop]]
            for kw in ['gene "|gene_name "', 'gene_id "', 'gene_biotype "']:
                temp = re.split(kw, features)
                if len(temp) > 1:
                    output.append(temp[1].split('"')[0])
                else:
                    output.append("unknown")
            if output[3] == "unknown":
                genes_skipped.add(output[4])
                continue
            if re.match(r'C..*orf\d', output[3]):
                genes_skipped.add(output[3])
                continue
            data.append(output)
    data = pd.DataFrame(data, columns = "chr start stop gene_name gene_id gene_biotype".split())
    log_message("chromosomes skipped: " + ", ".join(list(chrs_skipped)))
    log_message(f"{len(genes_skipped)} genes skipped due to name unknown or C..orf..")
    if args.GENES:
        genes_before = data.shape[0]
        biotypes = data["gene_biotype"].unique().tolist()
        if common := set(args.GENES) & set(biotypes):
            if missing := set(args.GENES) - set(biotypes):
                log_message("Biotypes not found in GTF:\n\t" + "\n\t".join(list(missing)))
            mask = data["gene_biotype"].isin(common)
            data = data[mask].copy()
        else:
            log_message("None of the biotypes are not found in the GTF.", exit_now=True)
        log_message(f"{data.shape[0]} of {genes_before} genes retained")
    if args.SORT:
        data["temp_chr"] = data["chr"].str.replace("chr","")
        data["temp_isdigit"] = [1 if x.isdigit() else 0 for x in data["temp_chr"]]
        for col in ["start", "stop"]:
            data[col] = data[col].astype(int)
        mask = data["temp_isdigit"] == 1
        autosomes = data[mask].copy()
        autosomes["temp_chr"] = autosomes["temp_chr"].astype(int)
        chromosomes_with_all_the_power = data[~mask].copy()
        autosomes.sort_values(by=["temp_chr", "start", "stop"], inplace=True)
        chromosomes_with_all_the_power.sort_values(by=["temp_chr", "start", "stop"], inplace=True)
        log_message("autosomes: " + ", ".join(autosomes["chr"].unique().tolist()))
        log_message("other chromosomes: " + ", ".join(chromosomes_with_all_the_power["chr"].unique().tolist()))
        data = pd.concat([autosomes, chromosomes_with_all_the_power], ignore_index=True, sort=False)
        data.drop(columns=["temp_chr", "temp_isdigit"], inplace=True)
    data.to_csv(outputfile, sep="\t", index=False)
    log_message(f"created {outputfile}")

def find_star_index(index=None):
    # takes basename of index, looks in indexstores
    index = os.path.basename(index)
    species = index.split(".")[0]
    if not species in constants.known_species:
        log_message(f"Unknown species {species}", exit_now=True)
    for dir in constants.index_stores:
        indexpath = os.path.join(dir, index)
        if os.path.exists(indexpath):
            if dir != constants.star.base_dir:
                log_message(
                    f"nice cp -r {indexpath} {constants.star.base_dir} &"
                )  # , exit_now=True)
            return indexpath, species
    temp = "\n".join(constants.index_stores)
    log_message(f"{index} not found in:\n{temp}", exit_now=True)

# STAR

def star_list_idx(args):
    print("\n".join(find_star_indexes().keys()))

def find_star_indexes():
    # return an orderedDict with index as key and dir as value
    index_path = OrderedDict()
    for dir in constants.index_stores:
        for index in map(
            os.path.basename, map(os.path.dirname, glob.glob(f"{dir}/*/SAindex"))
        ):
            if not index in index_path.keys():
                index_path[index] = dir
                if not index.split(".")[0] in constants.known_species:
                    log_message(
                        f"Caution: {dir}/{index} is not in constants.known_species"
                    )
        return index_path




if __name__ == "__main__":
    constants = define_constants()
    constants.star_indexes = find_star_indexes()
    args = define_args(constants)
    args.handle_args()
    fn = eval(args.func)
    fn(args.args)

