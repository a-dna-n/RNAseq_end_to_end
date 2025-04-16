#!/usr/bin/env python

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *

test_libraries(["cyclopts", "metapub"], exit_on_error = True)

import re
import pandas as pd
import cyclopts
from dataclasses import dataclass, KW_ONLY
from types import SimpleNamespace
from typing import Literal

cyc_app = cyclopts.App(help = "Utilities to manage workflow, studies, files.")
cyc_group = cyclopts.Group.create_ordered("workflow")

"""

def register():
    
    geo = "GSE237287"
    folder = "Polyak_24"

    check_dir_write_access(args.outputdir)

    srp = _pySRAdb_convert_ID(ID = geo, fn = "gse-to-srp")
    get_GEO_metadata(SimpleNamespace(study = geo, outputdir = f"{folder}/metadata",  parse=True, outputfileprefix = ""))
    pySRAdb_get_metadata(SimpleNamespace(study = srp, outputdir = f"{folder}/metadata",  parse=True, outputfileprefix = ""))
    get_NCBI_counts(SimpleNamespace(study = geo, outputdir =  f"{folder}/ncbi_recount"))
    geo_get_authors_files(input = [f"{folder}/metadata"], study = geo, outputdir = f"{folder}/author_files_in_GEO")


@dataclass
class study():
    _ : KW_ONLY
    key: str
    #base_dir: os.PathLike
    dirs = SimpleNamespace(
        exe = "${HOME}/bin/SPAdes-3.15.5/bin/rnaspades.py",
        tempdir = {"t3": "/media/2tb2/spades_temp"}.get(hostname, ""),
        mem = 60, #{"L490": 42, "t3": 60, "p43s": 20}.get(hostname, 0),
        sortmem = 4000 # {"L490": 4000, "t3": 4000, "p43s": 3000}.get(hostname, 0)
    )

class study_dirs:

    def __init__(self):
        hostname = execute_command("hostname")[0]
        read_type = SimpleNamespace(
            paired = "paired",
            single = "single"
        )
        read_types = [read_type.paired, read_type.single]
        # Literal conflicts with argparse
        # rnaseq.py genomeref: error: argument --species/-s: invalid choice: 'mouse' (choose from typing.Unpack[typing.Literal['human', 'mouse']])

        cpus = int(execute_command("grep -c ^processor /proc/cpuinfo")[0])
        known_species = ["human", "mouse"] #, "monkey", "vero"]
        _star_base_dir = "${HOME}/star"
        ref_dir = os.path.join(_star_base_dir, "ref")
        sortmem = 4000 #{"L490": 4000, "t3": 4000, "p43s": 3000}.get(hostname, 0)
        star = SimpleNamespace(
            base_dir = _star_base_dir,
            dummy_fastq = f"{_star_base_dir}/dummy.fastq",
            bin_dir = os.path.join(_star_base_dir, "scripts"),
            sortmem = int(sortmem / 2),
            sortcpus = int(cpus / 2),
            indexmem = 48000000000
        )
        coords_source_file = {
            str(species): os.path.join(ref_dir, f"gene_coords_{species}.txt")
            for species in known_species
        }
        index_stores = [_star_base_dir]
        index_stores += [
            x for x in ["/media/2tb2", "/mnt/p2"] if os.path.exists(x)
        ]  # keep in order

        featurecounts = SimpleNamespace(
            options = "--donotsort --fraction -M -O --extraAttributes gene_biotype,gene_name",
            constant_columns = "Geneid Chr Start End Strand Length gene_biotype gene_name".split(),
            unique_IDs = {
                species: f"{ref_dir}/ens.to.gene.if.unique.{species}" for species in known_species
            }
        )
        exit_if_startup_errors = False
        counts_sample_ID_file = "counts_sample_ids.txt"
        default_total_counts = "total_counts.txt"
        default_gene_metadata = "counts_gene_metadata.txt"
        default_merged_counts = "counts_combined.txt"
        default_counts_by_gene_name = "counts_by_gene_name.txt"
        geo_fetch_author_supp = "get_author_supp_files_from_GEO"
        geo_fetch_NCBI_counts = "get_NCBI_counts_human"
"""


@cyclopts.Parameter(name="*")
@dataclass
class _study_args:

    _: KW_ONLY
    
    study_key: str
    "Lastauthor_YY or Lastauthor_YY_MM"

    #parse: bool = True
    #"simplify/reformat the metadata"

@cyc_app.command(group=cyc_group)
def register(args: _study_args):
    """Make folder structure for a study
    """
    pass
    """
    #geo = 
    srp = _pySRAdb_convert_ID(ID = geo, fn = "gse-to-srp")
    get_GEO_metadata(SimpleNamespace(study = geo, outputdir = f"{folder}/metadata",  parse=True, outputfileprefix = ""))
    pySRAdb_get_metadata(SimpleNamespace(study = srp, outputdir = f"{folder}/metadata",  parse=True, outputfileprefix = ""))
    get_NCBI_counts(SimpleNamespace(study = geo, outputdir =  f"{folder}/ncbi_recount"))
    geo_get_authors_files(input = [f"{folder}/metadata"], study = geo, outputdir = f"{folder}/author_files_in_GEO")
    """

@cyc_app.command(group=cyc_group)
def html_to_md(*, inputfile: FileName, outputformat: Literal["gfm", "gfm-raw_html", "markdown", "markdown-raw_html"]="markdown-raw_html", outputfile: FileName | None = None):
    verify_that_paths_exist(inputfile)
    if not outputfile:
        outputfile = ".".join([inputfile, outputformat, "md"])
    exit_if_files_exist(outputfile)
    test_executables("pandoc", exit_on_error = True)
    execute_command(f"pandoc {inputfile} -o {outputfile} --to {outputformat} --from html --wrap=none") # --extract-media=figures")


if __name__ == "__main__":
    cyc_app()
