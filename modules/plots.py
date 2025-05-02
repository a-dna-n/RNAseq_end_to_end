#!/usr/bin/env python

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *

test_libraries(["cyclopts",], exit_on_error = True)

from types import SimpleNamespace
import re
import pandas as pd
import cyclopts
from dataclasses import dataclass, KW_ONLY
from typing import Annotated

cyc_app = cyclopts.App(help = "Functions for plots.") #, default_parameter=cyclopts.Parameter(consume_multiple=True))
cyc_group = cyclopts.Group.create_ordered("plots")
refdir = os.path.join(os.getenv("HOME"), "ptm/studies/ref")
ref_rna = {
    "human" : os.path.join(refdir, "trackplot_genome/Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz"),
    "mouse" : os.path.join(refdir, "trackplot_genome/Mus_musculus.GRCm39.109.sorted.gtf.gz")
}
ref_dna = {
    "human" : os.path.join(refdir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
    "mouse" : os.path.join(refdir, "Mus_musculus.GRCm39.dna.primary_assembly.fa")
}

def determine_species_from_ensembl_IDs(inputfile):
    verify_that_paths_exist(inputfile)
    with open(inputfile, "r") as tempin:
        line = next(tempin)
    if "ENSMUS" in line:
        return "mouse"
    if "ENSG" in line:
        return "human"
    log_message(f"Species not found in {inputfile}", fatal=True)

@cyclopts.Parameter(name="*")
@dataclass
class _junction_variant_list_args:

    _: KW_ONLY

    inputfile: FileName
    "Chr-start-stop coords + name."
 
    outputfile: FileName | None = None
    "Defaults to input + bed."

    overwrite: bool = False
    "Overwrite output."


def variant_def(x):
    chr, start, stop, gene = x["chr"], x["start"], x["stop"], x["gene"]

    return f"uniquely_mapped=1;multi_mapped=0;gene={gene};viewport={chr}:{start}-{stop};variant_name={gene}"


@cyc_app.command(group=cyc_group)
def output_junction_variant_list(args: _junction_variant_list_args):
    """Use gene coords for each study."""
    verify_that_paths_exist(args.inputfile)
    outputfile  = args.outputfile or args.inputfile + ".bed"
    if not args.overwrite:
         exit_if_files_exist(outputfile)
    """    
    with open(args.inputfile, "r") as tempin:
        with (open(regions, "w") as tempout:
            for line in tempin:
                coords = line.rstrip().split("\t")
                print(f"{coords[0]}:{coords[1]}-{coords[2]}", file=tempout)
    """

    data = pd.read_csv(args.inputfile, sep="\t", header=None)
    if len(data.columns) < 4:
        log_message(f"Invalid input", fatal=True)
    data = data[list(data.columns[:4])]
    data.columns = ["chr", "start", "stop",  "gene"]
    data["gene"] = data.apply(lambda x: variant_def(x), axis=1)
    data["ct"] = 1
    data["strand"] = "+"
    data.to_csv(outputfile, sep="\t", index=False, header = None)


@cyclopts.Parameter(name="*")
@dataclass
class _gtf_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    overwrite: bool = False
    "Overwrite output."


@cyc_app.command(group=cyc_group)
def make_local_gtf(args: _gtf_args):
    """Use gene coords for each study."""
    verify_that_paths_exist(args.inputfile)
    test_executables(["tabix", "bgzip"], exit_on_error = True)
    species = determine_species_from_ensembl_IDs(args.inputfile)
    gtf_out = args.inputfile + ".gtf"
    gz_out = gtf_out  + ".gz"
    if not args.overwrite:
        exit_if_files_exist([gtf_out, gz_out])
    cmd = f"tabix --regions {args.inputfile} {ref_rna[species]} > {gtf_out}" # 2> {gtf_out}.err")
    execute_command(cmd, output_command_first=True)
    cmd  = f"bgzip -k {gtf_out}" # 2> {gz_out}.err")
    execute_command(cmd, output_command_first=True)
    cmd  = f"tabix -p gff {gz_out}" # 2> {args.inputfile}.tabix.err")
    execute_command(cmd, output_command_first=True)


@cyclopts.Parameter(name="*")
@dataclass
class _gtf_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    overwrite: bool = False
    "Overwrite output."


@cyc_app.command(group=cyc_group)
def make_local_fa(args: _gtf_args):
    """Use gene coords for each study."""
    verify_that_paths_exist(args.inputfile)
    test_executables(["samtools"], exit_on_error = True)
    species = determine_species_from_ensembl_IDs(args.inputfile)
    regions = f"{args.inputfile}.regions"
    with open(args.inputfile, "r") as tempin:
        with open(regions, "w") as tempout:
            for line in tempin:
                coords = line.rstrip().split("\t")
                print(f"chr{coords[0]}:1-{int(coords[2])+10000}", file=tempout)

    fa_out = args.inputfile + ".fa"
    #gz_out = gtf_out  + ".gz"
    #if not args.overwrite:
    #    exit_if_files_exist([gtf_out, gz_out])
    cmd = f"samtools faidx --output {fa_out} --region-file {regions} {ref_dna[species]}" # 2> {gtf_out}.err")
    execute_command(cmd, output_command_first=True)
    if os.path.exists(fa_out):
        log_message(f"{fa_out} created")
    else:
        log_message(f"{cmd} failed")
        sys.exit()
    """
    cmd  = f"samtools faidx {fa_out}"
    execute_command(cmd, output_command_first=True)
    idx = fa_out + ".fai"
    if os.path.exists(idx):
        log_message(f"{idx} created")
    else:
        log_message(f"{cmd} failed")
        sys.exit()
    """

@cyclopts.Parameter(name="*")
@dataclass
class _junction_args:

    _: KW_ONLY

    inputfiles: Annotated[list[FileName], cyclopts.Parameter(consume_multiple=True)]

    "BAM input file(s)."

    overwrite: bool = False
    "Overwrite output."

@cyc_app.command(group=cyc_group)
def get_junctions(args: _junction_args):
    """Get junction counts from BAM files via regtools.
    
    output of regtools:

        15	98458558	98918721	JUNC00000001	1	?	98458558	98918721	255,0,0	2	74,49	0,460114
        15	98636357	98835733	JUNC00000009	1	?	98636357	98835733	255,0,0	2	88,19	0,199357
        15	98636357	98835811	JUNC00000010	1	?	98636357	98835811	255,0,0	2	88,19	0,199435

    Format for IGV:

        chr7	55019365	55142285	uniquely_mapped=180;multi_mapped=0;gene=EGFR	180	+
        chr7	55142437	55143304	uniquely_mapped=171;multi_mapped=0;gene=EGFR	171	+
        chr7	55143488	55146605	uniquely_mapped=228;multi_mapped=0;gene=EGFR	228	+

    """
    verify_that_paths_exist(args.inputfiles)
    test_executables(["regtools", "bgzip", "tabix"], exit_on_error = True)
    for inputfile in args.inputfiles:
        if not inputfile.endswith(".bam"):
            log_message("This function operates on BAM files.")
            continue
        junctions_out = inputfile + "_junctions"
        bed_out = junctions_out + ".bed"
        if not args.overwrite:
            exit_if_files_exist([junctions_out, bed_out])
        execute_command(f"regtools junctions extract {inputfile} -o {junctions_out} -M 25000 -s RFl", output_command_first=True)
        data = pd.read_csv(junctions_out, sep="\t", header=None)
        data.columns = "chr start stop desc count strand d1 d2 d3 d4 counts d5".split(" ")
        data["strand"] = "+"
        data["count"] = [sum(map(int,x.split(","))) for  x in data["counts"]]
        data["desc"] = [f"uniquely_mapped={ct};multi_mapped=0;" for ct in data["count"]]
        data["chr start stop desc count strand".split(" ")].to_csv(bed_out, sep="\t", header = None, index=False)

        execute_command(f"bgzip -k {bed_out}", output_command_first=True)
        execute_command(f"tabix -p bed {bed_out}.gz", output_command_first=True)


if __name__ == "__main__":
    cyc_app()


