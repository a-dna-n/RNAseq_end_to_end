#!/usr/bin/env python

"""
This module provides functions to retrieve RNA-seq data and metadata, launch alignments, counts, reassembly, etc.

It is not remotely intended to compete with Nextflow etc.

For some studies, there is data in SRA or GEO but not both, or only in ENA.
The reads may be in different place, but the ideal source is ENA since it lets us use aspera for download.

GEOparse and pySRAdb have excellent features, but neither knows about NCBI human recounts, and neither lets us choose which supplementary files to download depending on their contents. GEOparse only uses soft-formatted files, which are problematic in studies with human and mouse samples.

"""
"""
GSE154783, GSE228268  = c. sabaeus
"""
import sys
import os
from modules.tools import *

test_executables("free gunzip hostname pysradb samtools STAR wget".split(), exit_on_error = False)
test_libraries("numpy pandas pysradb".split(), exit_on_error = True) # GEOparse 

import re
import glob
from textwrap import dedent
from collections import OrderedDict, Counter
from types import SimpleNamespace
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
#import GEOparse
from modules.arg_def import *
from typing import Literal, Optional # TypeAlias, Union
#from lxml import html
#import requests
constants = define_constants()

from cyclopts import App as cyc_App, Group as cyc_Group, Parameter
from dataclasses import dataclass, KW_ONLY

cyc_app = cyc_App(help = "Functions for RNA-seq data and metadata.")
# function groups
metadata = cyc_Group.create_ordered("metadata")
buildrefs = cyc_Group.create_ordered("buildrefs")
align = cyc_Group.create_ordered("align")
counts = cyc_Group.create_ordered("counts")
assembly = cyc_Group.create_ordered("assembly")
bamfiles = cyc_Group.create_ordered("bamfiles")
utils = cyc_Group.create_ordered("utils")


test_cases_GEO = {
    # mouse - first study has multiple data types - should be filtered
    "mouse": "GSE275562 GSE263778 GSE286314".split(),
    # human - some have both  NCBI and author data
    "human": "GSE213519 GSE173475 GSE233661".split(),
    # chicken, zebrafish etc. - should not produce outputs beyond matrix files
    "xeno": "GSE278071 GSE87528 GSE283071 GSE276850".split(),
    # invalid ID but gets a response from GEO
    "invalid" : ["GSE239mmn"]
}

'''
# additional tests
PRJNA931290
SRP294329
SRP326996
SRP452987
GSE179462
SRP326996
'''



@cyc_app.command(group=metadata)
def _pySRAdb_convert_ID(*, ID: str, fn: str) ->str:
    """Simplistically calls pysradb to convert an ID using a specific function.
    
        > pysradb gsm-to-gse GSM4679562
        study_alias	study_accession
        GSE154783	SRP272683

    Args:
        ID (str): the ID to convert
        fn (str): the function to use

    Returns:
        str: ID in the specified namespace, e.g. GSE154783 in the example above
    """
    cmd = f"pysradb {fn} {ID}"
    log_message(cmd)
    temp = execute_command(cmd)
    prefix = fn[-3:].upper()
    if len(temp) == 2:
        new_IDs = temp[-1].rstrip().split()
        # new_id = temp[-1].rstrip().split()[1]
        for new_ID in new_IDs:
            if new_ID.startswith(prefix):
                log_message(f"{ID} => {new_ID}")
                return new_ID
    log_message(f"Unexpected output from pysradb:", *temp, sep="\n\t")
    return ""


@Parameter(name="*")
@dataclass
class _metadata_args:

    ID: str
    "SRA study ID, e.g. SRP294329"

    _: KW_ONLY
    
    outputdir: File_or_Dir
    "Folder where metadata files will be stored."

    outputfileprefix: str = constants.file_prefix.pysradb
    "Prefix for output file (base)names, to coerce nominal documentation of sources."

# def pySRAdb_get_metadata(ID: str, *, outputdir: File_or_Dir, outputfileprefix: str = constants.file_prefix.pysradb) -> pd.DataFrame:

@cyc_app.command(group=metadata)
def pySRAdb_get_metadata(args: _metadata_args) -> pd.DataFrame:
    """Find metadata and files available from SRA, GEO and and ENA.
    
    Args:
        args (_metadata_args): work in progress

    Returns:
        pd.DataFrame: if the ID is valid, returns a table populated with sample and fastq file info
      
    """    
    log_message(f"Retrieving metadata from SRA.")
    if not args.ID.startswith("SRP"):
        log_message(f"{args.ID} may not work. pySRAdb needs an SRP args.ID, e.g. SRP247494")

    # Get metadata for a study, save it as is.
    check_dir_write_access(args.outputdir)
    if os.path.basename(args.outputfileprefix) != args.outputfileprefix:
        log_message("Output file prefix cannot include directory", exit_now=True)

    if not constants.file_prefix.pysradb.lower() in args.outputfileprefix.lower():
        args.outputfileprefix = "_".join([args.outputfileprefix, constants.file_prefix.pysradb]) #{args.outputfileprefix}_pySRAdb"
    if not args.ID in args.outputfileprefix:
        args.outputfileprefix =  "_".join([args.outputfileprefix, args.ID])
    
    outputfile = os.path.join(args.outputdir, args.outputfileprefix + ".txt")
    if os.path.exists(outputfile):
        log_message(f"reading pySRAdb metadata from local file {outputfile}")
        return pd.read_csv(outputfile, sep="\t", header = 0, dtype=object), outputfile
    else:
        db = SRAweb()
        response = db.sra_metadata(args.ID, detailed=True)
        if response.shape[0]:
            response.to_csv(outputfile, sep="\t",  index=False)
            log_message(f"pySRAdb metadata for {args.ID} written to {outputfile}")
            return response, outputfile
        else:
            log_message(f"No metadata for {args.ID} via pySRAdb.")
            return None,  None


def parse_SRA_metadata(*, data: pd.DataFrame, ID: str, species: str|list[str], expt_type: str|list[str], outputfileprefix: str = constants.file_prefix.pysradb, overwrite: bool=False):
    """
    This code expects columns specific to RNA-seq.

    Subset metadata by species and assay, parse output separately.
    Only calls the metadata function fof pySRAdb, which can be used directly on the command line.
    if species is specified, returns a dict of dataframe(s) with species as the key
    otherwise returns a single dataframe
    """
    if isinstance(species, str):
        species = [species]
    if isinstance(expt_type, str):
         expt_type = [expt_type]
    if not  ID in outputfileprefix:
        outputfileprefix = "_".join([outputfileprefix, ID])
    outputfile = {(sp, ex): "_".join([outputfileprefix, constants.species_unalias[sp], ex])+".txt" for sp in species for ex in expt_type}
    if not overwrite:
        exit_if_files_exist(outputfile.values())
    #organism = {"human" : "Homo sapiens", "mouse" : "Mus musculus"}
    #if unk := set(species) - set(organism.keys()):
    #    log_message(f"Unknown species", *unk, sep="\n\t", exit_now=True)

    expt_type_columns = list(set(["strategy", "library_strategy"]) & set(data.columns))
    if len(expt_type_columns) == 0:
        log_message(f"No strategy or library_strategy column", exit_now=True)
    if len(expt_type_columns) > 1:
        log_message(f"Multiple strategy columns", exit_now=True)
    if "strategy" in data.columns:
        data.rename(columns = {"strategy": "library_strategy"}, inplace=True)

    columns_to_keep = ["run_total_bases", "run_total_spots", "organism_name", "library_layout", "study_accession"]

    verify_columns_present_in_dataframe(data=data, columns = columns_to_keep, source="SRA")
    
    # output summaries while filtering

    data_types = data["library_strategy"].unique().tolist()
    log_message(f"Types of data in this study:", *data_types, sep="\t")
    N_initial = data.shape[0]
    mask = data["library_strategy"].isin(expt_type)
    if data.shape[0] == 0:
        log_message(f"No samples match the experiment type.", exit_now=True)
    log_message(f"{data.shape[0]} rows out of {N_initial} were selected by type of data")

    species_in_study = data["organism_name"].unique().tolist()
    log_message(f"Species with data in this study:", *species_in_study, sep="\t")

    mask = data["organism_name"].isin(species)
    data = data[mask].copy()
    if data.shape[0] == 0:
        log_message(f"No samples have the selected species.", exit_now=True)
    log_message(f"{data.shape[0]} rows out of {N_initial} were selected by species")

    columns_to_drop = []

    def _check_experiment_alias(data):
        columns_to_drop = []
        if "experiment_alias" in data.columns:  # GSM
            mask = data["experiment_alias"].isna()
            tempdata = data[mask]
            if tempdata.shape[0]:
                log_message(f"{tempdata.shape[0]} missing run accessions:", *tempdata["run_accession"].unique().tolist(), sep=", ")
            tempdata = data[~mask]
            temp = list(tempdata["experiment_alias"])
            if len(temp) != len(set(temp)):
                log_message("Non-unique experiment_alias:", *set(temp), sep="\n\t")
            if "run_alias" in data.columns and list(data["run_alias"]) == list(
                    [f"{x}_r1" for x in data["experiment_alias"]]
                ):
                    columns_to_drop.append("run_alias")
        if columns_to_drop:
            return data.drop(columns = columns_to_drop)
        return data

    def SRA_simplify_exp_desc(x):
        temp = str(x).split(": ")
        if len(temp) == 1: return temp[0]
        return temp[1].split(";")[0]


    def _check_experiment_title(data):
        columns_to_drop = []
        if "experiment_title" in data.columns:
            for col in set(["experiment_desc", "library_name"]) & set(data.columns):
                if list(data[col]) == list(data["experiment_title"]):
                    columns_to_drop.append(col)
            data["experiment_title"] = data.apply(
                lambda x: SRA_simplify_exp_desc(x), axis=1
            )
        if columns_to_drop:
            return data.drop(columns = columns_to_drop)
        return data

    def _drop_empty_columns(data):
        if columns_to_drop := [x for x in data.columns if all(data[x].isna()) or all(data[x] == "missing")]:
            return data.drop(columns = columns_to_drop)
        return data

    def drop_unneeded(data):
        columns_to_drop = {col for col in data.columns if col[0:3] in ["aws", "gcp", "pub", "ncb"]}
        columns_to_drop  |= set("experiment_accession accession library_source library_selection instrument instrument_model instrument_model_desc organism_taxid total_size".split())
        if columns_to_drop := columns_to_drop & set(data.columns):
            return data.drop(columns = list(columns_to_drop))
        return data
   
    def delete_if_single_values(data, columns_to_keep):
        columns_to_drop = [col for col in data.columns if len(set(data[col]))==1]
        columns_to_drop = list(set(columns_to_drop) - set(columns_to_keep))
        if columns_to_drop := [x for x in data.columns if all(data[x].isna()) or all(data[x] == "missing")]:
            return data.drop(columns = columns_to_drop)
        return data
        
    def edit_ENA_fastq_URLs(data):
        if ena_columns := [col for col in data.columns if col.startswith("ena_fastq_")]:
            h_to_f = {}
            for col in ena_columns:
                if col.startswith("ena_fastq_ftp"):
                    data[col] = data[col].str.replace("era-fasp@fasp.sra.ebi.ac.uk:", "")
                elif col.startswith("ena_fastq_http"):
                    data[col] = data[col].str.replace("http://ftp.sra.ebi.ac.uk/", "")
                    h_to_f[col] =  col.replace("http", "ftp")
            if redundant := {h: f for h, f in h_to_f.items() if f in data.columns and data[h].tolist() == data[f].tolist()}:
                for h, f in redundant.items():
                    log_message(f"{h} is identical to {f}")
                return data.drop(columns = list(redundant.keys()))
        return data


    def edit_study_title(data):
        if "study_title" in data.columns:
            data["study_title"] = data["study_title"].str.replace(
                r"\s*\[RNA-Seq\]\s*", "", regex=True
            )
        return data

 
    data = _check_experiment_alias(data)
    data = _check_experiment_title(data)
    data = _drop_empty_columns(data)
    data = drop_unneeded(data)
    data = delete_if_single_values(data, columns_to_keep)
    data = edit_ENA_fastq_URLs(data)
    data = edit_study_title(data)

    data.insert(2, "read_length", 0)
    data.insert(3, "read_type", "")
    data.insert(4, "species", "")
    data["species"] = data["organism_name"].replace(constants.species_unalias)
    
    
    # calculate read length (inconsistent and variable)
    for col in ["run_total_bases", "run_total_spots"]:
        data[col] = data[col].astype(int)
    data["run_total_spots"] = [max(x,1) for x in data["run_total_spots"]]
    data["read_length"] = data["run_total_bases"] / data["run_total_spots"]
    data["read_length"]  = data["read_length"].round(0).astype(int)
    
    data.fillna("", inplace=True)
    data["library_layout"] = data["library_layout"].str.lower()

    if "ena_fastq_ftp_1" in data.columns and "ena_fastq_ftp_2" in data.columns:
        mask = (data["library_layout"] == constants.read_type.paired) | (
            (data["ena_fastq_ftp_1"].str.endswith(".gz"))
            & (data["ena_fastq_ftp_2"].str.endswith(".gz"))
        )
    else:
        mask = data["library_layout"] == constants.read_type.paired
    data.loc[mask, "read_length"] = data.loc[mask, "read_length"] / 2
    data.loc[mask, "read_type"] = constants.read_type.paired
    data.loc[~mask, "read_type"] = constants.read_type.single
    data["read_length"] = data["read_length"].round(0).astype(int)
    if constants.read_type.paired in set(data["read_type"]):
        data["read_length"] = data["read_length"].astype(str)
        mask = data["read_type"] == constants.read_type.paired
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
            mask = data["read_type"] == constants.read_type.paired
            for i in [1, 2]:
                data.loc[mask, f"sra_fastq_{i}"] = (
                    data["run_accession"] + f"_{i}.fastq"
                )
        else:
            mask = (data["library_layout"] == constants.read_type.paired) & (
                (data["ena_fastq_ftp_1"] == "") | (data["ena_fastq_ftp_2"] == "")
            )
            if mask.any():
                log_message("\nWarning: paired reads but single FTP file\n")
                data["sra_fastq_1"] = ""
                data["sra_fastq_2"] = ""
                mask = (data["library_layout"] == constants.read_type.paired) & (
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

    data = dedup_cols(data=data)
    
    output_files = []
    for sp in species:
            for exp in expt_type:
                mask = (data["organism_name"] == sp) & (data["library_strategy"] == exp)
                data_subset = data[mask]
                output_files.append(outputfile[(sp,exp)])
                if data_subset.shape[0]:
                    log_message(f"Writing {data_subset.shape[0]} rows with {exp} data for {sp} to {outputfile[(sp,exp)]}")
                    data_subset.to_csv(outputfile[(sp,exp)], sep="\t", index=False)
                    output_files.append(outputfile[(sp,exp)])
                else:
                    log_message(f"No metadata from SRA for {exp} data from {sp} in this study.")

    return output_files

'''

ENA can provide some sample info:

    read_experiment for accession PRJNA680934
        experiment_accession	run_accession	description	study_accession
        SRX9590782	SRR13150536	Illumina NovaSeq 6000 sequencing: GSM4946391: NC rep2 Homo sapiens RNA-Seq	PRJNA680934
        SRX9590784	SRR13150538	Illumina NovaSeq 6000 sequencing: GSM4946393: si-EIF5A2 rep1 Homo sapiens RNA-Seq	PRJNA680934
        SRX9590785	SRR13150539	Illumina NovaSeq 6000 sequencing: GSM4946394: si-EIF5A2 rep2 Homo sapiens RNA-Seq	PRJNA680934
        SRX9590783	SRR13150537	Illumina NovaSeq 6000 sequencing: GSM4946392: NC rep3 Homo sapiens RNA-Seq	PRJNA680934
        SRX9590786	SRR13150540	Illumina NovaSeq 6000 sequencing: GSM4946395: si-EIF5A2 rep3 Homo sapiens RNA-Seq	PRJNA680934
        SRX9590781	SRR13150535	Illumina NovaSeq 6000 sequencing: GSM4946390: NC rep1 Homo sapiens RNA-Seq	PRJNA680934
        
study:
    https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294329&result=study

    study_accession	description	secondary_study_accession
    PRJNA680934	RNA-Seq analysis EIF5A2 knockdown effect on ovarian cancer cells	SRP294329

taxon
    https://www.ebi.ac.uk/ena/portal/api/filereport?accession=9606&result=taxon
    tax_id	description
    9606	Homo sapiens

'''

def get_ENA_fastq_list(ID: str, *, outputdir: File_or_Dir, outputfileprefix: str = constants.file_prefix.ena): #outputfile: File_or_Dir | None = None) -> pd.DataFrame: #delete_temp_output: bool = False):
    
    # returns IDs  or  
    # gets list of ENA fastqs for ID
    # returns dataframe

    if not (ID.startswith("SRP") or ID.startswith("PRJ")):
        log_message(f"{ID} may not work. This ENA service needs an SRP or PRJ ID")
        log_message("For example, try SRP247494.")
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ID}&result=read_run&fields=run_accession,fastq_ftp"

    #https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP247494&result=read_run&fields=run_accession,fastq_ftp
    if  not constants.file_prefix.ena.lower() in outputfileprefix.lower():
        outputfileprefix = f"{outputfileprefix}_{constants.file_prefix.ena}"
    outputfile = f"{outputfileprefix}_{ID}.txt"
    outputfile = os.path.join(outputdir, outputfile)

    check_dir_write_access(outputdir)

    if os.path.exists(outputfile):
        if os.stat(outputfile).st_size ==  0:
            log_message(f"Error - {outputfile} is empty.",  exit_now=True)    
        log_message(f"Reading {outputfile} from previous query.")
        return pd.read_csv(outputfile, sep="\t", header  = 0), outputfile
    query = f"wget -q -O - '{url}'"
    if outputfile:
        query += f" | tee {outputfile}"
    log_message(f"Querying ENA for fastq files with:\n\t{query}")
    data = execute_command(query)
    if len(data):
        log_message(f"Retrieved {len(data)-1} records for fastq files")
        if outputfile:
            if os.stat(outputfile).st_size ==  0:
                log_message(f"Error - {outputfile} is empty.",  exit_now=True)    
            elif os.path.exists(outputfile):
                log_message(f"Check {outputfile}")
            else:
                log_message(f"{outputfile} missing", exit=True)
        data = [line.rstrip().split("\t") for line in data]
        data = pd.DataFrame(data[1:], columns =  data[0])
        return data, outputfile
    else:
        log_message(f"No results for {ID} with {query}")
        return pd.DataFrame(),  None

def reformat_ENA_fastq_list(data: pd.DataFrame, *, outputfile: File_or_Dir|None=None) -> pd.DataFrame:
    """
    Reformat so it matches SRA/pysradb columns. We also abbreviate the URL so it's easier to read when we choose  samples.

    run_accession	fastq_ftp
    SRR12879057	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/057/SRR12879057/SRR12879057_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/057/SRR12879057/SRR12879057_2.fastq.gz
    SRR12879059	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/059/SRR12879059/SRR12879059.fastq.gz
    SRR12879062	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/062/SRR12879062/SRR12879062_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/062/SRR12879062/SRR12879062_2.fastq.gz
    SRR12879068	ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/068/SRR12879068/SRR12879068.fastq.gz
    """

    #if not isinstance(data, pd.DataFrame):
    #    log_message(f"data is a {type(data)} instead of a dataframe")
    #    return None

    if outputfile:
        exit_if_files_exist(outputfile)
    columns_expected = ["run_accession", "fastq_ftp"]
    verify_columns_present_in_dataframe(data=data, columns = columns_expected, source="ENA fastq list")
    if unknown_columns := set(data.columns) - set(columns_expected):
        log_message("Dropping unknown columns:", *unknown_columns, sep="\n\t")
        data = data[columns_expected]
    data.rename(columns = {"fastq_ftp": "ena_fastq_ftp"}, inplace=True)
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
    data.reset_index(inplace=True)
    if outputfile:
        data.to_csv(outputfile, sep="\t", index=False)
        log_message(f"Reformatted ENA fastq list {outputfile}")
    return data


def test_get_ENA_fastq_list():
    data, outputfile = get_ENA_fastq_list(ID = "PRJNA931290", outputdir =  "test/ENA/PRJNA931290")
    data, outputfile = reformat_ENA_fastq_list(data = data, outputfile = outputfile)


def transpose_nested_list(*, data: list):
    # M x N list becomes N x M
    return [list(x) for x in zip(*data)]


def  test_URL(url):
    cmd = f"wget -o - --spider {url}"
    if response := execute_command(cmd):
        for line in response:
            if "Remote file exists" in line:
                return response
            if line.startswith("File ") and line.endswith(" exists."):
                return response
    return False

def get_GEO_metadata(ID: str, outputdir: File_or_Dir, list_files: bool = False)-> File_or_Dir: #, outputfileprefix: str = constants.file_prefix.geo

    if not re.match(r"GSE\d+$", ID):
        log_message("Species a series  ID e.g. GSE154891", exit_now=True)
    # input = GEO series ID, e.g. GSE154891
    # retrieves info for that GEO series:
    # GEO series home page - has links to author-submitted and NCBI counts (direct), indirect to matrix page
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154891

    # direct:https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162198/matrix/
    # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239889/matrix/

    check_dir_write_access(outputdir)
    #if os.path.basename(outputfileprefix) != outputfileprefix:
    #    log_message("Output file prefix cannot include directory", exit_now=True)
    series_ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{ID[:-3]}nnn/{ID}"
    matrix_page = f"{series_ftp_base}/matrix"
    matrix_file_urls = []
    for line in execute_command(f"wget -O  - -q {matrix_page}"):
        if "series_matrix.txt.gz" in line:
            temp = line.rstrip().split('"')[1]
            if temp.startswith("GSE") and temp.endswith(".gz"):
                 matrix_file_urls.append(f"{matrix_page}/{temp}")
            else:
                 log_message(f"Parsing error in {line} from {matrix_page}.", exit_now=True)
            # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239889/matrix/GSE239889-GPL24247_series_matrix.txt.gz
    local_files = []
    for url in matrix_file_urls:
        local_gz = os.path.join(outputdir, os.path.basename(url))
        local_file = local_gz.rstrip(".gz")
        if os.path.exists(local_file):
            pass
            #log_message(f"{local_file} present")
            #local_files.append(local_file)
        elif os.path.exists(local_gz):
            log_message(f"{local_gz} present - unip")
            execute_command(f"gunzip {local_gz}")
        else:
            
            #log_message(f"calling {cmd}")
            execute_command(f"wget -P {outputdir} {url}")
            execute_command(f"gunzip {local_gz}")
        if os.path.exists(local_file):
            local_files.append(local_file)
        else:
            log_message(f"Download or gunzip error for {url} - {local_gz} not found", exit_now=True)
    
    # retrieve filelist.txt if it exists
    file_list_url = f"{series_ftp_base}/suppl/filelist.txt"
    local_file_list = os.path.join(outputdir, os.path.basename(file_list_url))
    if os.path.exists(local_file_list):
        local_files.append(local_file_list)
    else:
        if test_URL(file_list_url):
            execute_command(f"wget --no-verbose --no-clobber -P {outputdir} {file_list_url}")
            local_files.append(local_file_list)
        else:
            log_message(f"{file_list_url} not found for {ID}")
    if list_files:
        log_message(f"metadata files for {ID}:", *local_files, sep="\n\t")
    return local_files

def get_NCBI_counts(GEO_ID: str, outputdir: File_or_Dir): #, destdir: File_or_Dir | None = None):
    """
    This function outputs a Bash script to retrieve files from NCBI. It does not retrieve the files.
    The Bash script doesn't retrieve the files directly. It outputs the commands so that they can be filtered on the fly.
    The commands are all generated by a function, so that this behavior is easy to change.
   
    GEO download pages have a consistent format. It's silly to parse them repeatedly.
    
    Example:
        Source: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162198

        /geo/download/?type=rnaseq_counts&acc=GSE162198&format=file&file=GSE162198_raw_counts_GRCh38.p13_NCBI.tsv.gz
        /geo/download/?type=rnaseq_counts&acc=GSE162198&format=file&file=GSE162198_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz
        /geo/download/?type=rnaseq_counts&acc=GSE162198&format=file&file=GSE162198_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
        /geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz

    This breaks down to:

        base_URL="https://www.ncbi.nlm.nih.gov/geo/download"
        {base_URL} / ? type=rnaseq_counts & acc={GEO_ID} & format=file & file= {GEO_ID}_raw_counts_GRCh38.p13_NCBI.tsv.gz
        {base_URL} / ? type=rnaseq_counts & acc={GEO_ID} & format=file & file= {GEO_ID}_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz
        {base_URL} / ? type=rnaseq_counts & acc={GEO_ID} & format=file & file= {GEO_ID}_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
        {base_URL} / ? type=rnaseq_counts                & format=file & file= Human.GRCh38.p13.annot.tsv.gz
    """
    script = os.path.join(outputdir, f"{constants.geo_fetch_NCBI_counts}_{GEO_ID}.sh")
    if os.path.exists(script):
        log_message(f"{script} exists for {GEO_ID}")
        return
    
    check_dir_write_access(outputdir)
    # may not exist
    
    download = "https://www.ncbi.nlm.nih.gov/geo/download/?"
    acc = download + f"acc={GEO_ID}"
    response = execute_command(f'wget -O - "{acc}"', splitlines=False)
    if not "Gene Expression Omnibus" in response:
        log_message(f"Failed to locate {GEO_ID} with {acc}")
        return
    if not "NCBI-generated data" in response:
        log_message(f"No NCBI-generated data located for {GEO_ID} with {acc}.")
        return
    #if not test_URL(acc):
    #   invalid test - GEO responds even with GSE342nnn
    #    log_message(f"Download page not found for {GEO_ID} at {acc}", exit_now=True)
    output = [dedent(f"""
        #!/bin/bash
        # set -e
        # 
        # This script fetches gene expression counts and derived values (TPM and FPKM) generated by GEO staff, assuming they exist for the series (GSE..)
        # Here is an explanation of the effort: https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html
        # The comments below were copied from the page linked above.
        # You can find an example here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162198
        # 
        #my_dir=$(dirname $(readlink -f ${{BASH_SOURCE}}))
        my_dir=$(dirname ${{BASH_SOURCE}})
        base_URL=
        GEO_series="{GEO_ID}"
        acc="${{GEO_series}}&"
        download_file="https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file"
        download_series_file="${{download_file}}&acc=${{GEO_series}}"
        
        go() {{
            file=$1
            url=$2
            echo ""
            echo "# $file"
            echo ""
            local_file=${{my_dir}}/${{file}}

            if [[ -e ${{local_file}} ]]; then 
                echo "${{local_file}} exists" >&2
            else
                if [[ -e ${{local_file}}.gz ]]; then
                    echo "${{local_file}}.gz exists" >&2
                else
                    echo "wget --no-clobber --quiet -O ${{local_file}}.gz \\"${{url}}&file=${{file}}.gz\\" && sleep 1"
                fi
                echo "gunzip --keep ${{local_file}}.gz"
            fi
        }}

        # Series RNA-seq raw counts matrix
        # Raw count matrices may be suitable for input into differential expression analysis tools.
        go "${{GEO_series}}_raw_counts_GRCh38.p13_NCBI.tsv" ${{download_series_file}}

        # Series RNA-seq normalized counts matrix
        # Normalized count matrices may be suitable for analyzing and visualizing gene expression abundance.
        go "${{GEO_series}}_norm_counts_FPKM_GRCh38.p13_NCBI.tsv" ${{download_series_file}}
        go "${{GEO_series}}_norm_counts_TPM_GRCh38.p13_NCBI.tsv" ${{download_series_file}}

        # Human gene annotation table
        go Human.GRCh38.p13.annot.tsv ${{download_file}}
    """).lstrip()]
    with open(script, "w") as tempout:
        print(*output, sep="\n", file=tempout)
        os.chmod(script, 0o755)
        log_message(f"Created {script}")

    """
    baseurl = f"{download}/?acc={GEO_ID}"

    ncbi_counts = []
    #response = execute_command(f'wget -q -O - "{baseurl}"', splitlines=False)
    response = requests.get(baseurl)
    """
    """
    >>> type(response.content)
    <class 'bytes'>
    >>> type(response.text)
    <class 'str'>
    >>> 
    """
    """
    tree = html.fromstring(response.content)
    urls = [x[2] for x in list(tree.iterlinks())]
    ncbi_counts = list(filter(lambda row: "rnaseq_counts" in row, urls))
    ncbi_counts = [x.replace("/geo/download/","") for x in ncbi_counts]
    
    for line in :
        if m := re.search(r'href="(.*?NCBI.tsv.gz)"', line):
            # /geo/download/?type=rnaseq_counts&amp;acc=GSE215024&amp;format=file&amp;file=GSE215024_raw_counts_GRCh38.p13_NCBI.tsv.gz
            ncbi_counts.append(m.groups()[0])
    if not ncbi_counts:
        log_message(f"No counts files were found for {GEO_ID} at {baseurl}")
        return
    """
    """
    for f in ncbi_counts:
        file = f.split("=")[-1]
        output.append(f"wget --no-clobber -q -O {file} 'https://www.ncbi.nlm.nih.gov{f}'")
    """

def matrix_to_dataframes(inputfile):
    all_data  = []
    entities = set()
    with open(inputfile, "r") as tempin:
        for line in tempin:
            tabs = line[1:].split("\t", 1)
            if len(tabs) > 1:
                entity, info = tabs[0].split("_", 1)
                entities.add(entity)
                all_data.append([entity,  info, tabs[1].rstrip().replace('"', '').split("\t")])
    data_by_entity = {}
    for ent in entities:
        data_by_entity[ent] = list(filter(lambda row: row[0] == ent, all_data))

    sample_data = data_by_entity["Sample"]
    test = [x[1] for x in sample_data]
    exclude = set()
    skip_if_match = "channel_count data_row_count platform_id status type growth_protocol".split(" ")
    skip_if_contain = "contact extract_protocol growth_protocol _date".split(" ")
    for i in test:
        if any(i in s for s in skip_if_match):
            exclude.add(i)
            continue
        if any(s in i for s in skip_if_contain):
            exclude.add(i)
            continue
    sample_data  = [x for x in sample_data if not x[1] in exclude]

    for i, row in enumerate(sample_data):
        sample_data[i][1] = sample_data[i][1].replace("_ch1", "")
    sample_data = [x[1:] for x in sample_data]
    variable_columns = []
    common_columns = []
    for row in sample_data:
        values = Counter(*row[1:])
        if len(values.keys()) == 1:
            val = list(values.keys())[0]
            if not (val in ["None", "NONE", "0"]):
                common_columns.append(row)
        else:
            variable_columns.append(row)

    info = variable_columns + common_columns
    for i, row in enumerate(sample_data):
        if row[0].startswith("characteristics"):
            h = set([x.split(": ")[0] for x in row[1]])
            if len(h) == 1:
                sample_data[i][0] = h.pop().replace(" ","_")
                for j, r in enumerate(row[1]):
                    sample_data[i][1][j] = r.split(": ")[1]
    sample_data = unique_row_headers(sample_data)
    sample_data = [[x[0]] + x[1] for x in sample_data]
    temp = transpose_nested_list(data=sample_data)
    sample_metadata = pd.DataFrame(temp[1:], columns=temp[0])

    series_data = [x[1:] for x in data_by_entity["Series"]]
    row_headers = [x[0] for x in series_data]
    exclude = set()
    skip_if_match = "status contributor".split(" ")
    skip_if_contain = "contact _date platform_".split(" ")
    for i in row_headers:
        if any(i in s for s in skip_if_match):
            exclude.add(i)
            continue
        if any(s in i for s in skip_if_contain):
            exclude.add(i)
            continue
    series_data  = [[x[0]] + x[1] for x in series_data if not x[0] in exclude]
    id_mapping = {"geo_accession" : "study", "SRA" : "SRP", "BioProject"  : "bioproject",  "pubmed_id" : "Pubmed"}
    ids_found = {}

    for i, x in enumerate(series_data):
        if mapping := id_mapping.get(x[0], None):
            ids_found[mapping] = x[1]
        match x[0]:
            case "sample_taxid":
                if org := constants.species_tax_code.get(x[1], ""):
                    series_data.append(["species", org])
                    if not "organism" in sample_metadata.columns:
                        log_message("Rescue org/species here")
                    #series_data[i] = ["species", org]
            case "relation":
                temp = x[1].split(": ", 1)
                if len(temp) > 1:
                    match temp[0]:
                        case "SRA":
                            temp[1] = temp[1].split("=")[1]
                        case "BioProject":
                            temp[1] = temp[1].split("/")[-1]
                    series_data.append(temp)
                    if mapping := id_mapping.get(temp[0], None):
                        ids_found[mapping] = temp[1]
    for key in reversed(id_mapping.values()):
        if key in ids_found:
            sample_metadata.insert(0, key, ids_found[key])

    series_metadata = pd.DataFrame(series_data, columns=["field", "value"])

    return sample_metadata, series_metadata

    
def parse_GEO_matrix_files(inputfiles, *, species: str|list[str]|None=None, expt_type: str|list[str]|None=None, outputfileprefix: str = constants.file_prefix.geo, overwrite: bool=False):
    # input = name of a series matrix file
    # output = a dataframe with sample-specific info and a dataframe with common entries
    # to do: check species & exp type
    #return_data =[]
    return_files = []
    if species is not None and isinstance(species, str):
        species = [species]
    if expt_type is not None and isinstance(expt_type, str):
         expt_type = [expt_type]
    for matrix_file in inputfiles:
        if not "matrix.txt" in matrix_file:
            log_message(f"Skipping {matrix_file}")
            continue
        matrix_file = matrix_file.replace(".gz", "")
        gz = f"{matrix_file}.gz"
        if (not os.path.exists(matrix_file)) and os.path.exists(gz):
            execute_command(f"gunzip {gz}")
        verify_that_paths_exist(matrix_file)
        log_message(f"\nParsing {matrix_file}")
        sample_metadata, series_metadata = matrix_to_dataframes(matrix_file)

        #sample_metadata.to_csv("sample_metadata", sep="\t", index=False)
        #series_metadata.to_csv("series_metadata", sep="\t")
        #sys.exit()

        verify_columns_present_in_dataframe(data=sample_metadata, columns=["organism", "library_strategy"], source=matrix_file)
        assays_in_data = sample_metadata["library_strategy"].unique().tolist()
        log_message(f"Types of data in this series or subseries:", *assays_in_data, sep="\t")
        if expt_type:
            N_initial = sample_metadata.shape[0]
            mask = sample_metadata["library_strategy"].isin(expt_type)
            sample_metadata = sample_metadata[mask].copy()
            if sample_metadata.shape[0] == 0:
                log_message(f"No samples match the specified experiment type(s).")
            else:
                log_message(f"{sample_metadata.shape[0]} rows out of {N_initial} were selected by type of data")
        organisms_in_data = sample_metadata["organism"].unique().tolist()
        log_message(f"Organisms found in this series or subseries:", *organisms_in_data, sep="\t")
        if len(organisms_in_data) > 1:
             log_message(f"Warning: multiple species were found in the same matrix file {matrix_file}.")
        if species:
            N_initial = sample_metadata.shape[0]
            mask = sample_metadata["organism"].isin(species)
            sample_metadata = sample_metadata[mask].copy()
            if sample_metadata.shape[0] == 0:
                log_message(f"No samples match the specified organism(s).")
            else:
                log_message(f"{sample_metadata.shape[0]} rows out of {N_initial} were selected by organism")
        for org in sample_metadata["organism"].unique().tolist():
            for assay in sample_metadata["library_strategy"].unique().tolist():
                outputfile = "_".join([matrix_file.replace(".txt",""), constants.species_unalias[org], assay]) + ".txt"
                mask = (sample_metadata["organism"] == org) & (sample_metadata["library_strategy"] == assay)
                tempdata = sample_metadata[mask].copy()
                if columns_to_drop := [col for col in tempdata.columns if set(tempdata[col]) in [{'NONE'}, {""}]]:
                    #vals = {col : set(tempdata[col]) for col in tempdata.columns}
                    #if columns_to_drop := [col for col, vals in vals.items() if len(vals)==1  and (vals[0] in ["NONE",  ""])]:
                    tempdata.drop(columns = columns_to_drop, inplace=True)

                for col in tempdata.columns:
                    splits = [str(x).split(": ") for x in tempdata[col]]
                    test = splits[0][0]
                    if all([(len(x) == 2 and x[0] == test) for x in splits]):
                        tempdata[col] = [x[1] for x in splits]
                        tempdata.rename(columns = {col : test}, inplace=True)

                tempdata.to_csv(outputfile, sep="\t", index=False)
                log_message(f"Matrix {matrix_file} converted to:\n\t{outputfile}\nfor {org} and {assay}\n")
                #return_data.append(tempdata)
                return_files.append(outputfile)
        # output parsed version of series matrix file, script to retrieve author suppplementary files
        outputfile = matrix_file.replace(".txt","_table.txt")
        series_metadata.to_csv(outputfile, sep="\t", index=False)
        #print(f"{matrix_file} {org} {assay} {tempdata.shape[0]}")
        log_message(f"Series matrix {matrix_file} parsed, reformatted, written to {outputfile}\n")
        #geo_fetch_author_supp = "get_author_supplementary_files_from_GEO_"
        #return_data.append(tempdata)
        return_files.append(outputfile)

    return return_files

    """
        all_data  = []
        entities = set()
        with open(matrix_file, "r") as tempin:
            for line in tempin:
                tabs = line[1:].split("\t", 1)
                if len(tabs) > 1:
                    entity, info = tabs[0].split("_", 1)
                    entities.add(entity)
                    all_data.append([entity,  info, tabs[1].rstrip().replace('"', '').split("\t")])
        data_by_entity = {}
        for ent in entities:
            data_by_entity[ent] = list(filter(lambda row: row[0] == ent, all_data))
        sample_data = data_by_entity["Sample"]
        test = [x[1] for x in sample_data]
        exclude = set()
        skip_if_match = "channel_count data_row_count platform_id status type growth_protocol".split(" ")
        skip = "channel_count data_row_count platform_id status type growth_protocol".split(" ")
        skip_if_contain = "contact extract_protocol growth_protocol _date".split(" ")
        for i in test:
            if any(i in s for s in skip_if_match):
                exclude.add(i)
                continue
            if any(s in i for s in skip_if_contain):
                exclude.add(i)
                continue
        sample_data  = [x for x in sample_data if not x[1] in exclude]
 
        for i, row in enumerate(sample_data):
            sample_data[i][1] = sample_data[i][1].replace("_ch1", "")
        sample_data = [x[1:] for x in sample_data]
        variable_columns = []
        common_columns = []
        for row in sample_data:
            values = Counter(*row[1:])
            if len(values.keys()) == 1:
                val = list(values.keys())[0]
                if not (val in ["None", "NONE", "0"]):
                    common_columns.append(row)
            else:
                variable_columns.append(row)

        info = variable_columns + common_columns
        for i, row in enumerate(sample_data):
            if row[0].startswith("characteristics"):
                h = set([x.split(": ")[0] for x in row[1]])
                if len(h) == 1:
                    sample_data[i][0] = h.pop().replace(" ","_")
                    for j, r in enumerate(row[1]):
                        sample_data[i][1][j] = r.split(": ")[1]
        sample_data = unique_row_headers(sample_data)
        sample_data = [[x[0]] + x[1] for x in sample_data]
        temp = transpose_nested_list(data=sample_data)
        out = pd.DataFrame(temp[1:], columns=temp[0])
        out.to_csv("yo", sep="\t")
        print(b)
        sys.exit()
        
    return  pd.DataFrame(b[1:], columns=b[0])
    """
"""

        script = f"{constants.geo_fetch_author_supp}_{ID}.sh"
        geo_find_supplementary_files(files=matrix_files, outputfile=script)

"""

def geo_find_supplementary_files(
    *, input: File_or_Dir|list[File_or_Dir], outputdir: File_or_Dir, ID: str, overwrite: bool = False
):
    # input = list of series matrix files
    # output = Bash script to retrieve any supplementary files, or None if None
    #verify_that_paths_exist(files)
    check_dir_write_access(outputdir)
    files = []
    if isinstance(input, File_or_Dir):
        input = [input]
    for i in input:
        if os.path.isdir(i):
            if temp := glob.glob(f"{i}/*series_matrix.txt"):
                files.extend(temp)
            else:
                log_message(f"No matrix files found in {i}")
        elif os.path.isfile(i):
            if i.endswith("series_matrix.txt"):
                files.append(i)
            else:
                log_message(f"{i} is an invalid name for a matrix file")
        else:
            log_message(f"Skipping {i}")
    if not files:
        log_message(f"No matrix files")
        return
    script = f"{constants.geo_fetch_author_supp}_{ID}.sh"
    script = os.path.join(outputdir, script)
    if not overwrite:
        exit_if_files_exist(script)

    output = []
    output.append(dedent(f"""
        #!/bin/bash
        set -e
        # fetch supplementary files posted by authors of {ID} in GEO
        # destination dir is where this Bash script resides
        dest_dir=$(dirname ${{BASH_SOURCE}})

        get() {{
            url=$1
            local_file=${{dest_dir}}/$(basename ${{url}})
            echo "wget -nc -nv -O ${{local_file}} ${{url}}"
        }}

    """))
        
    ftp_files_found = []
    # log_message("Ignore warnings like: Execution of wget -nc --spider ftp://... failed.")
    for matrix_file in files:
        if matrix_file.endswith(".gz"):
            matrix_file = matrix_file.replace(".gz", "")
            if not os.path.exists(matrix_file):
                execute_command(f"gunzip {matrix_file}.gz")
        tempoutput = []
        with open(matrix_file, "r") as tempin:
            lines = tempin.readlines()
        for line in lines:
            if line.startswith("!Series_supplementary_file"):
                url = line.split('"')[1]
                if url in ftp_files_found:
                    continue
                if url.startswith("ftp"):
                    if response := test_URL(url):
                        b = os.path.basename(url)
                        for g in response:
                            if f"SIZE {b}" in g:
                                temp = g.rstrip().split()[-1]
                                tempoutput.append(f"# {b} = {int(temp):,} bytes")
                                break
                        tempoutput.append(f"get {url}\n")
                        ftp_files_found.append(url)
                else:
                    log_message(f"Malformed FTP URL in {line} from {matrix_file}", exit_now=True)
        if tempoutput:
            #output.append(f"# source: {matrix_file}\n")
            output += tempoutput
    if ftp_files_found:
        with open(script, "w") as tempout:
            print(*output, sep="\n", file=tempout)
        os.chmod(script, 0o755)
        log_message(f"created {script}")
        return script
    else:
        log_message(f"No FTP file URLs found.")
        return None

def geo(args):

    convert_functions = set("""
        gse-to-gsm gse-to-srp gsm-to-gse gsm-to-srp gsm-to-srr gsm-to-srs gsm-to-srx 
        srp-to-gse srp-to-srr srp-to-srs srp-to-srx srr-to-gsm srr-to-srp srr-to-srs 
        srr-to-srx srs-to-gsm srs-to-srx srx-to-gsm srx-to-srp srx-to-srr srx-to-srs
    """.split())

    direct_to_srp = {x.split("-to-")[0] : x for x in convert_functions if x.endswith("srp")}
    direct_to_gse = {x.split("-to-")[0] : x for x in convert_functions if x.endswith("gse")}
    skip = [x for x in convert_functions if x.startswith("gse")]
    convert_functions -= set(skip)
    convert_functions -= set([x for x in convert_functions if  x.split("-to-")[0] in direct_to_gse])
    intermediates = set([x for x in convert_functions if  x.split("-to-")[1] in direct_to_gse])
    indirects = set([x.split("-to-")[0] for x in intermediates])
    paths_to_try = {x: [i for i in intermediates  if i.split("-to-")[0] == x] for x in indirects}

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
        get_NCBI_counts(GEO_ID=ID)
        outputfile = f"{ID}{outputfileext}"
        if os.path.exists(outputfile):
            log_message(f"skipping {ID} - output file {outputfile} exists")
            continue
        matrix_files = get_GEO_series_metadata_files(geo_id=ID)

        script = f"{constants.geo_fetch_author_supp}_{ID}.sh"
        geo_find_supplementary_files(files=matrix_files, outputfile=script)
        if problematic := [x for x in matrix_files if not ("_series_matrix" in x)]:
            #temp = "\n".join(problematic)
            log_message(
                "_series_matrix not found for matrix files:", *problematic, sep="\n\t", exit_now=True
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
                SRP = srp_ids_found[0]
                log_message(f"Caution - multiple SRP IDs were found for {id_succession[0]}:",  *srp_ids_found,  sep="\n\t")
                log_message(f"Using {SRP}")

        for i, file in enumerate(matrix_files):

            sample_metadata, series_metadata = matrix_to_dataframes(file)
            sample_metadata.insert(0, "study", id)
            sample_metadata.insert(1, "SRP", SRP)
            sample_metadata.insert(2, "source_", subseries[i])
            sample_metadata = dedup_cols(data=sample_metadata)
            sample_metadata.to_csv(output_files[i], sep="\t", index=False)
            log_message(f"Metadata written to {output_files[i]}")
        if SRP:
            SRP = _pySRAdb_convert_ID(fn="gse-to-srp", id=id)
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
                        cmd = f"{__file__} join -i {outputfile} {sra_output} -c geo_accession {sra_column} -o {temp}"
                        execute_command(cmd)
                        if os.path.exists(temp):
                            log_message(f"GEO + SRA => {temp}")
                        else:
                            log_message(f"Failed to create {temp} with:\n{cmd}")
                else:
                    log_message(f"No matched columns in {outputfile} and {sra_output}")


def test_get_GEO_series_metadata_files():

    for species, studies in test_cases_GEO.items():
        base_output_dir = f"test/GEO_metadata/{species}"
        for ID in studies:
            files = get_GEO_metadata(ID = ID, outputdir = f"{base_output_dir}/{ID}", list_files=True)
            files = parse_GEO_matrix_files(files, species = ["Homo sapiens", "Mus musculus"], expt_type = "RNA-Seq")

def test_get_NCBI_counts():

    for species, studies in test_cases_GEO.items():
        base_output_dir = f"test/NCBI_counts/{species}"
        for ID in studies:
            get_NCBI_counts(GEO_ID = ID, outputdir = f"{base_output_dir}/{ID}")
          
def  test_geo_find_supplementary_files():

    for species, studies in test_cases_GEO.items():
        base_input_dir = f"test/GEO_metadata/{species}"
        base_output_dir = f"test/GEO_author_files/{species}"
        for ID in studies:
            geo_find_supplementary_files(input=f"{base_input_dir}/{ID}", outputdir=f"{base_output_dir}/{ID}", ID=ID, overwrite=True) 

def bam_species(args):
    if args.OUTPUTFILE:
        exit_if_files_exist(args.OUTPUTFILE)
    verify_that_paths_exist(args.BAMFILES)
    species = get_species_for_bamfiles(bamfiles=args.BAMFILES)
    with open(args.OUTPUTFILE, "w") if args.OUTPUTFILE else sys.stdout as tempout:
        for bam in args.BAMFILES:
            print(f"{bam}\t{species[bam]}", file=tempout)

def bash_header(*, flags: str = "set -e"):
    return dedent(f"""
        #!/bin/bash
        {flags}
        log() {{ echo "$*" >&2 ; }}
        die() {{ echo "$*" >&2; exit 1 ; }}
    """).lstrip()


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

r'''
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
        x: os.path.realpath(re.sub(r"\.bam$", "", x)) for x in bamfiles
    }
    for d in outputdir_by_sample.values(): check_dir_write_access(d)

    read_type = {x: get_read_type_from_a_bamfile(bamfile=x) for x in bamfiles}
    fastq_files = {}
    fastq_file = {}
    for bam in bamfiles:
        prefix = re.sub(r"\.bam$", "", bam)
        if read_type[bam] == constants.read_type.paired:
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
        if read_type[bam] == constants.read_type.paired:
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
'''

def handle_commands_and_output(*,
    commands: str | list[str], outputfile: File_or_Dir|None = None, single_line_command: str|None =  None, execute: bool = False
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
    pass
    r'''
    from tempfile import mkstemp

    # needs coords and either gene bams or non-gene BAMs, but not both
    test_executables("regtools")
    outputfileext = args.EXT
    bamfiles = args.BAMFILES
    verify_that_paths_exist(bamfiles)
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
    '''

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
    if read_type == constants.read_type.single:
        command += " --single"
    elif read_type != constants.read_type.paired:
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
                output.append(f"gunzip --keep -q {ref_file}.gz")
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
        log_message(f"""
            zcat {gtf}.gz > {gtf}
               or
            gunzip {gtf}
        """, exit_now=True)
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
        print(*output, sep="\n", file=tempout)

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
            print(*output, sep="\n", file=tempout)
        os.chmod(args.SCRIPT, 0o755)
        log_message(f"created {args.SCRIPT}")
    else:
        print(*output, sep="\n")

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


def bash_featurecounts(*, gtf: File_or_Dir, cpus: int, options: str):
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


def bash_fasterq(*, cpus: int, read_type: constants.read_types):
    bash_function_name = "fastqs_sra"
    opts = {constants.read_type.single: "", constants.read_type.paired: "--split-files"}.get(read_type)
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
            sampleopts="--outFile_or_DirPrefix $prefix."
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
            log_message(f"Header row found in {inputfile}:", *temp, sep="\n\t")
    data = pd.read_csv(inputfile, sep="\t", header=header, dtype=object)
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
        read_type = constants.read_type.single
    elif len(fastq_columns) == 2:
        read_type = constants.read_type.paired
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
        if read_type == constants.read_type.single:
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
        if read_type == constants.read_type.single:
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

def check_if_ref_files_match(*, dna: File_or_Dir, rna: File_or_Dir, species: constants.known_species):

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
    """
    if args.RNASPADESCOORDS:
        verify_that_paths_exist(args.RNASPADESCOORDS)
        test_executables(rnaspades.exe)
    """
    outputdir = args.OUTPUTDIR
    counts_dir = os.path.join(outputdir, "counts")
    check_dir_write_access([outputdir, counts_dir])
    if glob.glob(f"{counts_dir}/*"):
        log_message(f"\nCaution: files found in {counts_dir}\n", exit_now=True)
    bamdestdir = args.BAMDESTDIR or ""
    if bamdestdir:
        check_dir_write_access(bamdestdir)
    if args.ABRA2COORDS:
        abra2dir = os.path.join(outputdir, "abra2")
        check_dir_write_access(abra2dir)
    """
    if args.RNASPADESCOORDS:
        rnaspadesdir = os.path.join(outputdir, "rnaspades")
        check_dir_write_access(rnaspadesdir)
        # rnaspadesdir = os.path.join(outputdir, "rnaspades")
        # temp = os.path.basename(os.path.realpath(outputdir))
        # rnaspadesdir = check_dir_write_access(dir = os.path.join(rnaspadesdir, temp))
    """
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

    if bamdestdir and args.ABRA2COORDS: # or args.RNASPADESCOORDS):
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
    """
    if args.RNASPADESCOORDS:
        srcdir = bamdestdir or os.path.realpath(outputdir)
        output.append(f"cp -s {srcdir}/*.bam {srcdir}/*.bai {rnaspadesdir}")
        output.append(f"mv {args.RNASPADESCOORDS} {rnaspadesdir}")
        output.append(f"cd {rnaspadesdir}")
        output.append(
            f"{__file__} rnaspades -b *.bam -r {os.path.basename(args.RNASPADESCOORDS)} -z --exec"
        )
        output.append(f"cd {os.path.realpath(outputdir)}")
    """
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
        #if outputfile and os.sep in outputfile:
        if outputfile and os.path.basename(args.outputfile) != args.outputfile:
            log_message(
                "outputfile cannot include dir if outputdir is specified", exit_now=True
            )
        if totals and os.path.basename(totals) != totals:
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


def rank_values(args):
    verify_that_paths_exist(args.INPUTFILE)
    if args.OUTPUTFILE and not args.OVERWRITE:
        exit_if_files_exist(args.OUTPUTFILE)

    temp = pd.read_csv(args.INPUTFILE, sep="\t", header=0, comment="#", nrows=1)
    if args.COLUMN:
        verify_columns_present_in_dataframe(
            data=temp, columns=args.COLUMN, source=args.INPUTFILE
        )
        if not args.NEWCOLUMN:
            newcolumn = f"rank_by_{args.COLUMN}"
    else:
        if args.SORT:
            log_message("There is no sense in sorting by the rank", exit_now=True)
        newcolumn = args.NEWCOLUMN or "order"
    verify_columns_absent_in_dataframe(
            data=temp,  columns=newcolumn, source=args.INPUTFILE
        )
    data = pd.read_csv(args.INPUTFILE, sep="\t", header=0, comment="#")

    if args.COLUMN:
        if args.SORT:
            data.sort_values(by=args.COLUMN, ascending = args.ORDER == "ascending", inplace=True)
        data[newcolumn] = data[args.COLUMN].rank(method="min").astype(int)
    else:
        data[newcolumn] = list(range(1, data.shape[0]+1))
    data.to_csv(args.OUTPUTFILE if args.OUTPUTFILE else sys.stdout, sep="\t", index=False)


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
        STAR --genomeDir {args.INDEX} --readFilesIn {constants.star.dummy_fastq} --runThreadN 4 --outFile_or_DirPrefix ${constants.star.dummy_fastq}. --genomeLoad Remove --runMode alignReads"
    """))
    # --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 25 --seedMultimapNmax 250 --seedPerReadNmax 250 --alignTranscriptsPerReadNmax 1000 --limitOutSJoneRead 100 --alignWindowsPerReadNmax 250 
    output = output.dedent()
    with open(outputscript, "w") if outputscript else sys.stdout as tempout:
        print(output, file=tempout)
    if outputscript:
        os.chmod(outputscript, 0o755)
        log_message(f"created {outputscript}")
        #print(outputscript)

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


def dedup_cols(data: pd.DataFrame):
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

r'''
def minimap_output_file(inputfile: File_or_Dir):  # =None, ext=None):
    outputfile = re.sub(r"\.gz$", "", inputfile)
    outputfile = re.sub(".fasta", "", outputfile)
    # outputfile = re.sub(r"\.fa$", "", outputfile)
    # if ext:
    #    outputfile = f'{outputfile}.{ext.lstrip(".")}'
    return outputfile


def bash_minimap(*, cpus: int, bases_loaded: str, genome: str, junctions: File_or_Dir):
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
'''

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
    # (inputfile: File_or_Dir, *, outputfile: File_or_Dir|None=None, sort: bool = False, subset=["protein_coding"]):
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

    # (inputfile: File_or_Dir, *, outputfile: File_or_Dir|None=None, sort: bool = False, subset=["protein_coding"]):
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
    print(*constants.star_indexes, sep="\n") #"\n".join(find_star_indexes().keys()))

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


def get_GEO_series_metadata_files(geo_id: str, *, unzip: bool=True)-> list[str]:
    if not re.match(r"GSE\d+$", geo_id):
        log_message("Species a series  ID e.g. GSE154891", exit_now=True)
    # input = GEO series ID, e.g. GSE154891
    # retrieves info for that GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154891
    # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239889/matrix/

    matrix_files = []
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}"
    log_message(f"fetching {url}")
    ftp = ""
    for line in execute_command(f"wget -O - {url}"):
        if m := re.search(r'.+(ftp:.*?)".+Series Matrix File', line):
            ftp = m.groups()[0]
            #print(ftp)
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
            # print(ftp)
            break
    if not ftp:
        log_message(f"No series matrix file found in {url}")
        return

    local_flat_files = []
    local_gz_files = []
    local_gz = os.path.basename(ftp)
    local_file = local_gz.replace(".gz", "")
    if os.path.exists(local_file):
        log_message(f"{local_file} exists")
        matrix_files.append(local_file)
    elif os.path.exists(local_gz):
        log_message(f"{local_gz} already retrieved")
        matrix_files.append(local_gz)
    else:
        cmd = f"wget --quiet --no-clobber {ftp} 2> /dev/null"
        log_message(f"executing {cmd}")
        execute_command(cmd)
        if os.path.exists(local_gz):
            log_message(f"Retrieved {local_gz}")
            matrix_files.append(local_gz)
        else:
            log_message(f"Failed to retrieve {local_gz} with {cmd}")
    if unzip:
        for i, file in enumerate(matrix_files):
            if file.endswith(".gz"):
                cmd = f"gunzip -q {file}"
                log_message(f"executing {cmd}")
                execute_command(cmd)
                file.replace(".gz", "")

    return matrix_files

def metadata(args):
    """
    Procedure:
        1. retrieve the available metadata from SRA, GEO and ENA, saving original files by study then source.
        2. filter. reformat and separate by species
        3. make scripts to retrieve supplementary files and NCBI recounts if available.
        4. merge sample details with fastq details, choose a subset of the samples, and make a simple manifest
        5. make a read-length-specific index if we don't have one already
        6 start the alignments
    """
    studydir = os.path.realpath(args.STUDYDIR)
    subdir = {s: os.path.join(studydir, s) for s in args.SOURCES}
    check_dir_write_access(studydir + list(subdir.values()))

    id_prefix_by_source = {"GEO": "GEO", "SRP": "SRA", "PRJ": "ENA"}
'''    
    if "SRA" in source_for_db = SRAweb()

    g2=GEOparse.get_GEO(geo="GSE162198", silent=True, destdir="workspace/GSE162198")
    response = db.sra_metadata(ID, detailed=True)
    for source in source_for_ids.keys():
        tempdir = subdir[source]
 
     match source:
        case "GEO":

       >>> g2.relations
        {'BioProject': ['https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680934'], 'SRA': ['https://www.ncbi.nlm.nih.gov/sra?term=SRP294329']}

        >>> gse.relations
        {'BioProject': ['https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1001347']}

            return "Bad request"
        case "SRA":
            #d = _pysradb(["SRP294329"], outputdir = "workspace/pysradb",  outputfileprefix ="sra_", species = ["human", "mouse"],
            #             expt_type = ["RNA-seq"])
            #data, file = pySRAdb_get_metadata(ID = "SRP294329", outputdir = "workspace/pysradb/SRP294329", outputfileprefix = "pySRAdb")
            #data = parse_SRA_metadata(data=data, ID = "SRP294329", species = ["Homo sapiens", "Mus musculus"], expt_type = "RNA-Seq", #outputfileprefix = file.replace(".txt","")) #os.path.join("workspace/pysradb/SRP294329", "pySRAdb"))
            
            data, file = pySRAdb_get_metadata(ID = "GSE239889", outputdir = "workspace/pysradb/SRP294329", outputfileprefix = "pySRAdb")
            data = parse_SRA_metadata(data=data, ID = "GSE239889", species = ["Homo sapiens", "Mus musculus"], expt_type = "RNA-Seq", outputfileprefix = file.replace(".txt","")) #os.path.join("workspace/pysradb/SRP294329", "pySRAdb"))
     case "ENA":
            # we need SRP or PRJ for  ENA
            data = get_ENA_fastq_list(ID = ID, outputfile)
            data = parse_ENA_fastq_metadata(data, outputfile)

            data = get_ENA_fastq_list(ID = "PRJNA931290", outputfile =  "workspace/pysradb/PRJNA931290")
            data = reformat_ENA_fastq_list(data = data, outputfile = "workspace/pysradb/PRJNA931290.txt")
            # if ena columns are not in SRA metadata,  get two columns from ENA fastq manifest.
        sys.exit()
   
         case _:
            log_message(f"unknown source {source}", exit_now=True)
  
'''


if __name__ == "__main__":

    test_geo_find_supplementary_files()
    test_get_NCBI_counts()
    test_get_GEO_series_metadata_files()
    test_get_ENA_fastq_list()
    sys.exit()
    constants.star_indexes = find_star_indexes()
    args = define_args(constants)
    args.handle_args()
    fn = eval(args.func)
    fn(args.args)
    #cyc_app()
