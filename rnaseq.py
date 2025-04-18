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
test_libraries("numpy pandas pysradb metapub cyclopts".split(), exit_on_error = True) # GEOparse 
import cyclopts
"""
from modules import metadata, expression, viz_bokeh

cyc_app = cyc_App(help = "Functions for cyno RNA-seq analysis.")
cyc_plots = cyc_Group.create_ordered("plots")

cyc_app.update(metadata.metadata_app)
cyc_app.update(expression.app)
#cyc_app.update(expression.expression_app)
cyc_app.update(viz_bokeh.app)


"""


import re
import glob
from textwrap import dedent
from collections import OrderedDict, Counter
from types import SimpleNamespace
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
#import GEOparse
#from modules.arg_def import *
from typing import Literal, Optional # TypeAlias, Union
from modules.constants import *
#from lxml import html
#import requests
constants = define_constants()

from cyclopts import App as cyc_App, Group as cyc_Group, Parameter
from dataclasses import dataclass, KW_ONLY, field

cyc_app = cyc_App(help = "Functions for RNA-seq data and metadata.")
# function groups
cyc_groups = {}
cyc_groups["data"] = cyc_Group.create_ordered("data")
"""
buildrefs = cyc_Group.create_ordered("buildrefs")
align = cyc_Group.create_ordered("align")
counts = cyc_Group.create_ordered("counts")
assembly = cyc_Group.create_ordered("assembly")
bamfiles = cyc_Group.create_ordered("bamfiles")
"""
cyc_groups["utils"] = cyc_Group.create_ordered("utils")


# def get_NCBI_counts(*, study: str, outputdir: File_or_Dir): #, destdir: File_or_Dir | None = None):

@Parameter(name="*")
@dataclass
class _data_args:

    _: KW_ONLY

    study: str
    "SRA or GEO study ID e.g. SRP294329, GSE117552"
    
    outputdir: File_or_Dir
    "Folder where metadata files will be stored."

    #outputfileprefix: str = constants.file_prefix.pysradb
    #"Prefix for output file (base)names, to coerce nominal documentation of sources."

# def get_NCBI_counts(*, study: str, outputdir: File_or_Dir): #, destdir: File_or_Dir | None = None):



@cyc_app.command(group=cyc_groups["data"])
def get_NCBI_counts(args: _data_args): #, destdir: File_or_Dir | None = None):
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
    script = os.path.join(args.outputdir, f"{constants.geo_fetch_NCBI_counts}_{args.study}.sh")
    if os.path.exists(script):
        log_message(f"{script} exists for {args.study}")
        return
    
    check_dir_write_access(args.outputdir)
    # may not exist
    
    download = "https://www.ncbi.nlm.nih.gov/geo/download/?"
    acc = download + f"acc={args.study}"
    response = execute_command(f'wget -O - "{acc}"', splitlines=False)
    if not "Gene Expression Omnibus" in response:
        log_message(f"Failed to locate {args.study} with {acc}")
        return
    if not "NCBI-generated data" in response:
        log_message(f"No NCBI-generated data located for {args.study} with {acc}.")
        return
    #if not test_URL(acc):
    #   invalid test - GEO responds even with GSE342nnn
    #    log_message(f"Download page not found for {args.study} at {acc}", exit_now=True)
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
        GEO_series="{args.study}"
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
    baseurl = f"{download}/?acc={args.study}"

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
        log_message(f"No counts files were found for {args.study} at {baseurl}")
        return
    """
    """
    for f in ncbi_counts:
        file = f.split("=")[-1]
        output.append(f"wget --no-clobber -q -O {file} 'https://www.ncbi.nlm.nih.gov{f}'")
    """


'''    
def metadata(args):
    """
    Procedure:
        1. retrieve the available metadata from SRA, GEO and ENA, saving original files by study then source.
        2. filter, reformat and separate by species
        3. make scripts to retrieve supplementary files and NCBI recounts if available.
        4. merge sample details with fastq details, choose a subset of the samples, and make a simple manifest
        5. make a read-length-specific index if we don't have one already
        6 start the alignments
    """
    studydir = os.path.realpath(args.STUDYDIR)
    subdir = {s: os.path.join(studydir, s) for s in args.SOURCES}
    check_dir_write_access(studydir + list(subdir.values()))

    id_prefix_by_source = {"GEO": "GEO", "SRP": "SRA", "PRJ": "ENA"}

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
    #constants.star_indexes = find_star_indexes()
    cyc_app()


