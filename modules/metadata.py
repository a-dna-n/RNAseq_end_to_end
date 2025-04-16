#!/usr/bin/env python

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *
from modules.constants import *
test_libraries(["cyclopts",], exit_on_error = True)

from types import SimpleNamespace
import re
import pandas as pd
import cyclopts
from dataclasses import dataclass, KW_ONLY, field
from typing import Literal

constants = define_constants()
cyc_app = cyclopts.App(help = "Functions for metadata.")
cyc_group = cyclopts.Group.create_ordered("metadata")

test_cases_GEO = {
    # mouse - first study has multiple data types - should be filtered
    "mouse": "GSE275562 GSE263778 GSE286314".split(),
    # human - some have both NCBI and author data
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

@cyc_app.command(group=cyc_group)
def test_get_NCBI_counts():

    for species, studies in test_cases_GEO.items():
        base_output_dir = f"test/NCBI_counts/{species}"
        for ID in studies:
            get_NCBI_counts(GEO_ID = ID, outputdir = f"{base_output_dir}/{ID}")

@cyc_app.command(group=cyc_group)          
def  test_geo_get_authors_files():

    for species, studies in test_cases_GEO.items():
        base_input_dir = f"test/GEO_metadata/{species}"
        base_output_dir = f"test/GEO_author_files/{species}"
        for ID in studies:
            geo_get_authors_files(input=f"{base_input_dir}/{ID}", outputdir=f"{base_output_dir}/{ID}", ID=ID, overwrite=True) 





@cyc_app.command(group=cyc_group)
def geo_get_authors_files(
    *, input: list[File_or_Dir], outputdir: File_or_Dir, study: str, overwrite: bool = False
):
    # input = list of series matrix files
    # output = Bash script to retrieve any supplementary files, or None if None
    #verify_that_paths_exist(files)
    check_dir_write_access(outputdir)
    files = []
    inputs = []
    for i in input:
        if os.path.isdir(i):
            print("isdir")
            inputs.extend(glob.glob(f"{i}/*series_matrix.txt"))
        elif "*" in i:
            print("wildcard")
            inputs.extend(glob.glob(i))
        else:
            print("other")
            inputs.append(i)
    for i in inputs:
        """
        if os.path.isdir(i):
            if temp := glob.glob(f"{i}/*series_matrix.txt"):
                files.extend(temp)
            else:
                log_message(f"No matrix files found in {i}")
        """
        if os.path.isfile(i):
            if i.endswith("series_matrix.txt"):
                files.append(i)
            else:
                log_message(f"{i} is an invalid name for a matrix file")
        else:
            log_message(f"Skipping {i}")
    if not files:
        log_message(f"No matrix files")
        return
    script = f"{constants.geo_fetch_author_supp}_{study}.sh"
    script = os.path.join(outputdir, script)
    if not overwrite:
        exit_if_files_exist(script)

    output = []
    output.append(dedent(f"""
        #!/bin/bash
        set -e
        # fetch supplementary files posted by authors of {study} in GEO
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





@cyc_app.command(group=cyc_group)
def get_GEO_series_metadata_files(geo_id: str, *, unzip: bool=True)-> list[str]:
    # "Get metadata from GEO"
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



@cyclopts.Parameter(name="*")
@dataclass
class _metadata_args:

    _: KW_ONLY

    study: str
    "SRA or GEO study ID e.g. SRP294329, GSE117552"

    outputdir: File_or_Dir
    "Folder where metadata files will be stored."

    outputfileprefix: str = constants.file_prefix.pysradb
    "Prefix for output file (base)names, to coerce nominal documentation of sources."

    parse: bool = True
    "simplify/reformat the metadata"

# def pySRAdb_get_metadata(ID: str, *, outputdir: File_or_Dir, outputfileprefix: str = constants.file_prefix.pysradb) -> pd.DataFrame:

@cyc_app.command(group=cyc_group)
def pySRAdb_get_metadata(args: _metadata_args) -> None: #pd.DataFrame:
    """Find metadata and files available from SRA, GEO and and ENA.
    
    Args:
        args (_metadata_args): work in progress

    Returns:
        Nothing - writes output to file(s)
    """    
    #  pd.DataFrame: if the ID is valid, returns a table populated with sample and fastq file info
    log_message(f"Retrieving metadata from SRA.")
    if args.study.startswith("GSE"):
        log_message("Attempt to obtain SRP from pysradb")
        args.study = _pySRAdb_convert_ID(ID = args.study, fn = "gse-to-srp")
    if not args.study.startswith("SRP"):
        log_message(f"{args.study} may not work. pySRAdb needs an SRP study, e.g. SRP247494")

    # Get metadata for a study, save it as is.
    check_dir_write_access(args.outputdir)
    if os.path.basename(args.outputfileprefix) != args.outputfileprefix:
        log_message("Output file prefix cannot include directory", exit_now=True)

    if not constants.file_prefix.pysradb.lower() in args.outputfileprefix.lower():
        args.outputfileprefix = "_".join([args.outputfileprefix, constants.file_prefix.pysradb]) #{args.outputfileprefix}_pySRAdb"
    if not args.study in args.outputfileprefix:
        args.outputfileprefix =  "_".join([args.outputfileprefix, args.study])
    
    outputfile = os.path.join(args.outputdir, args.outputfileprefix + ".txt")
    if os.path.exists(outputfile):
        log_message(f"reading pySRAdb metadata from local file {outputfile}")
        data = pd.read_csv(outputfile, sep="\t", header = 0, dtype=object)
    else:
        db = SRAweb()
        data = db.sra_metadata(args.study, detailed=True)
        if data is None or data.shape[0] == 0:
            log_message(f"No metadata for {args.study} via pySRAdb.")
            return # None,  None
        data.to_csv(outputfile, sep="\t",  index=False)
        log_message(f"pySRAdb metadata for {args.study} written to {outputfile}")
    if args.parse:
        parse_SRA_metadata(data=data, ID = args.study, species = ["Homo sapiens", "Mus musculus", "Chlorocebus sabaeus"], expt_type = "RNA-Seq", outputfileprefix = outputfile.replace(".txt",""), overwrite=True)

    #return data

    #return pd.read_csv(outputfile, sep="\t", header = 0, dtype=object), outputfile


@cyc_app.command(group=cyc_group)
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

    verify_columns_present(data=data, columns = columns_to_keep, source="SRA")
    
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
 
    data = _drop_empty_columns(data)
    data = _check_experiment_alias(data)
    data = _check_experiment_title(data)
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

    # return output_files

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

"""
@dataclass
class enafastqs:
help_text ="Get a ist of ENA fastqs for ID(s) like PRJNA627881"


    IDs: list[str], required=True
    "Identifier(s)."


    keeptempfiles: bool = False
    "Keep temp output files."


    ext: str = "ena.txt"
    "Output file extension"

"""


@cyc_app.command(group=cyc_group)
def get_ENA_fastq_list(args: _metadata_args) -> None: #ID: str, *, outputdir: File_or_Dir, outputfileprefix: str = constants.file_prefix.ena): #outputfile: File_or_Dir | None = None) -> pd.DataFrame: #delete_temp_output: bool = False):
    "Get the list of fastq files from ENA for a study."
    
    # returns IDs  or  
    # gets list of ENA fastqs for args.study
    # returns dataframe

    if not (args.study.startswith("SRP") or args.study.startswith("PRJ")):
        log_message(f"{args.study} may not work. This ENA service needs an SRP or PRJ args.study")
        log_message("For example, try SRP247494.")
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={args.study}&result=read_run&fields=run_accession,fastq_ftp"

    #https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP247494&result=read_run&fields=run_accession,fastq_ftp
    if  not constants.file_prefix.ena.lower() in args.outputfileprefix.lower():
        args.outputfileprefix = f"{args.outputfileprefix}_{constants.file_prefix.ena}"
    outputfile = f"{args.outputfileprefix}_{args.study}.txt"
    outputfile = os.path.join(args.outputdir, outputfile)

    check_dir_write_access(args.outputdir)

    if os.path.exists(outputfile):
        if os.stat(outputfile).st_size ==  0:
            log_message(f"Error - {outputfile} is empty.",  exit_now=True)    
        log_message(f"Reading {outputfile} from previous query.")
        return pd.read_csv(outputfile, sep="\t", header  = 0), outputfile
    query = f"wget -q -O - '{url}'" #query += f" | tee {outputfile}" jj- failed
    log_message(f"Querying ENA for fastq files with:\n\t{query}") 
    data = execute_command(query)
    if len(data) > 1:
        log_message(f"Retrieved {len(data)-1} records for fastq files")
        if outputfile:
            with open(outputfile, "w") as tempout:
                print(*data, sep="\n", file=tempout)
            log_message(f"ENA fastq details {args.study} written to {outputfile}")
        """
        if outputfile:
            if os.stat(outputfile).st_size ==  0:
                log_message(f"Error - {outputfile} is empty.",  exit_now=True)    
            elif os.path.exists(outputfile):
                log_message(f"Check {outputfile}")
            else:
                log_message(f"{outputfile} missing", exit=True)
        """
        """
        data = [line.rstrip().split("\t") for line in data]
        data = pd.DataFrame(data[1:], columns =  data[0])
        data.to_csv(outputfile, sep="\t", index=False)
        log_message(f"ENA fastq details {args.study} written to {outputfile}")
        #return data, outputfile
        """
    else:
        log_message(f"No results for {args.study} with {query}")
        return
        #return pd.DataFrame(),  None
    if args.parse:
        data = [line.rstrip().split("\t") for line in data]
        data = pd.DataFrame(data[1:], columns =  data[0])
        data = reformat_ENA_fastq_list(data, outputfile = outputfile.replace(".txt", ".parsed.txt"))



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
    verify_columns_present(data=data, columns = columns_expected, source="ENA fastq list")
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

@cyc_app.command(group=cyc_group)
def test_get_ENA_fastq_list():
    data, outputfile = get_ENA_fastq_list(ID = "PRJNA931290", outputdir =  "test/ENA/PRJNA931290")
    data, outputfile = reformat_ENA_fastq_list(data = data, outputfile = outputfile)


@cyc_app.command(group=cyc_group)
def get_GEO_metadata(args: _metadata_args) -> None: #pd.DataFrame:*, study: str, outputdir: File_or_Dir, list_files: bool = False)-> File_or_Dir: #, outputfileprefix: str = constants.file_prefix.geo
    "Get metadata from GEO a study."

    if not re.match(r"GSE\d+$", args.study):
        log_message("Species a series  args.study e.g. GSE154891", exit_now=True)
    # input = GEO series args.study, e.g. GSE154891
    # retrieves info for that GEO series:
    # GEO series home page - has links to author-submitted and NCBI counts (direct), indirect to matrix page
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154891

    # direct:https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162198/matrix/
    # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239889/matrix/

    check_dir_write_access(args.outputdir)
    #if os.path.basename(outputfileprefix) != outputfileprefix:
    #    log_message("Output file prefix cannot include directory", exit_now=True)
    series_ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{args.study[:-3]}nnn/{args.study}"
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
        local_gz = os.path.join(args.outputdir, os.path.basename(url))
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
            execute_command(f"wget -P {args.outputdir} {url}")
            execute_command(f"gunzip {local_gz}")
        if os.path.exists(local_file):
            local_files.append(local_file)
        else:
            log_message(f"Download or gunzip error for {url} - {local_gz} not found", exit_now=True)
    
    # retrieve filelist.txt if it exists
    file_list_url = f"{series_ftp_base}/suppl/filelist.txt"
    local_file_list = os.path.join(args.outputdir, args.study + "_" + os.path.basename(file_list_url))
    if os.path.exists(local_file_list):
        local_files.append(local_file_list)
    else:
        if test_URL(file_list_url):
            execute_command(f"wget --no-verbose --no-clobber -O {local_file_list} {file_list_url}")
            #execute_command(f"wget --no-verbose --no-clobber -P {args.outputdir} {file_list_url}")
            local_files.append(local_file_list)
        else:
            log_message(f"{file_list_url} not found for {args.study}")
    #if list_files:
    if args.parse:
        local_files += parse_GEO_matrix_files(local_files, species = ["Homo sapiens", "Mus musculus", "Chlorocebus sabaeus"], expt_type = "RNA-Seq")
    log_message(f"metadata files for {args.study}:", *local_files, sep="\n\t")
    #return local_files



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

    
@cyc_app.command(group=cyc_group)
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

        verify_columns_present(data=sample_metadata, columns=["organism", "library_strategy"], source=matrix_file)
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
        geo_get_authors_files(files=matrix_files, outputfile=script)

"""

@cyclopts.Parameter(name="*")
@dataclass
class geo_md_args:

    studies: list[str]
    "Identifier(s)."

    ext: str = "_geo.txt"
    "Output file suffix"


@cyc_app.command(group=cyc_group)
def geo(args: geo_md_args) -> None:
    """Get metadata from GEO"""

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
        
    outputfileext = args.ext
    # outputfile = {id : f"{id}.{outputfileext}" for id in ids}
    # exit_if_files_exist(list(outputfile.values()))
    for ID in args.studies:

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
        geo_get_authors_files(files=matrix_files, outputfile=script)
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

@cyc_app.command(group=cyc_group)
def test_get_GEO_series_metadata_files():

    for species, studies in test_cases_GEO.items():
        base_output_dir = f"test/GEO_metadata/{species}"
        for study in studies:
            files = get_GEO_metadata(study = study, outputdir = f"{base_output_dir}/{study}", list_files=True)
            files = parse_GEO_matrix_files(files, species = ["Homo sapiens", "Mus musculus", "Chlorocebus sabaeus"], expt_type = "RNA-Seq")

@cyc_app.command(group=cyc_group)
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


@cyclopts.Parameter(name="*")
@dataclass
class concat_file_args:

    inputfiles: list[FileName]
    "Input file(s)."


    outputfile: File_or_Dir | None = None
    "Output file. Stdout if omitted."

@cyc_app.command(group=cyc_group)
def concat_files(args: concat_file_args):
    "Concatenate files as dataframes, handle different columns"
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


        
@cyclopts.Parameter(name="*")
@dataclass
class _join_args:

    _: KW_ONLY
    
    inputs: list[str]
    "Two or more files."

    columns: list[str] = field(default_factory=lambda: ["Sample_geo_accession", "experiment_alias"])
    "Column(s) to use, as one (if identical) or two comma-separated lists."

    outputfileprefix: str = constants.file_prefix.pysradb
    "Prefix for output file (base)names, to coerce nominal documentation of sources."

    dedup: bool = False
    "Deduplicate columns after merging."

    outputfile: File_or_Dir | None = None
    "Output file"

    method: Literal[*["inner", "outer", "left", "right"]] = "inner"
    "How to join."


@cyc_app.command(group=cyc_group)
def join_files(args: _join_args) -> None: #pd.DataFrame:
    Nfiles = len(args.inputs)
    if Nfiles < 2:
        log_message("Specify multiple input files", exit_now=True)
    verify_that_paths_exist(args.inputs)
    outputfile = args.outputfile
    if outputfile:
        exit_if_files_exist(outputfile)
    method = args.method
    columns = args.columns
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
    for i, inputfile in enumerate(args.inputs):
        data = pd.read_csv(inputfile, sep="\t", header=0, dtype=object)
        verify_columns_present(
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
    if args.dedup:
        output = dedup_cols(data=output)
    output.to_csv(outputfile if outputfile else sys.stdout, sep="\t", index=False)

'''
@cyc_app.command(group=cyc_group)
def parse_fastq_manifest(*, i: FileName) -> None: #pd.DataFrame:
    """Parse fastq alignment manifest, output SRR and assigned label
    """
    inputfile = i
    verify_that_paths_exist(inputfile)
    with open(i, "r") as  tempin:
        for line in tempin:
            temp = line.rstrip().split("\t")
            if len(temp) < 2:
                log_message("Invalid format", fatal=True)
            files = [t for t in temp if "/" in t]
            samples = list(set(temp) - set(files))
            if len(files) ==  0 or len(samples) == 0:
                log_message("Invalid format", fatal=True)
            file = os.path.basename(files.pop())
            if file.startswith("SRR"):
                file = file.split("_")[0]
            else:
                log_message("Not a file from SRA - update this function")
            print(*[file, samples[0]],  sep="\t")
'''


@cyc_app.command(group=cyc_group)
def parse_fastq_manifest(*, i: FileName) -> None: #pd.DataFrame:
    """Parse fastq alignment manifest, output SRR and assigned label
    """
    inputfile = i
    verify_that_paths_exist(inputfile)
    with open(i, "r") as  tempin:
        for line in tempin:
            temp = line.rstrip().split("\t")
            if len(temp) < 2:
                log_message("Invalid format", fatal=True)
            files = [t for t in temp if "/" in t]
            samples = list(set(temp) - set(files))
            if len(files) ==  0 or len(samples) == 0:
                log_message("Invalid format", fatal=True)
            file = os.path.basename(files.pop())
            if file.startswith("SRR"):
                file = file.split("_")[0]
            else:
                log_message("Not a file from SRA - update this function")
            print(*[file, samples[0]],  sep="\t")




if __name__ == "__main__":
    cyc_app()


