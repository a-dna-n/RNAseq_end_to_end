#!/usr/bin/env python

# CAUTION - work in progress - arg handling is only partially migrated to cyclopts

import sys
import os
from types import SimpleNamespace
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *
from modules.constants import *
from collections import OrderedDict
import glob
#test_libraries(["cyclopts",], exit_on_error = True)

from types import SimpleNamespace
import re
import pandas as pd
import cyclopts
from dataclasses import dataclass, KW_ONLY
from typing import Literal

constants = define_constants()

cyc_app = cyclopts.App(help = "Functions for RNA-seq analyses from scratch.")
cyc_group_align = cyclopts.Group.create_ordered("pipeline")
cyc_group_refs = cyclopts.Group.create_ordered("Preprocessing.")


@cyc_app.command(group=cyc_group_refs)
def genomeref(args):
    if args.script:
        if "{species}" in args.script:
            args.script = args.script.replace("{species}", args.species)
        if not args.overwrite:
            exit_if_files_exist(args.script)
    dna = constants.ref.dna[args.species]
    rna = constants.ref.rna[args.species]
    output = []
    output.append(bash_header())
    output.append(dedent(f"""
        destdir={args.outputdir}
        dna={dna}
        rna={rna}
        wget --no-clobber -P $destdir $dna
        wget --no-clobber -P $destdir $rna
    """))
    if args.script:
        with open(args.script, "w") as tempout:
            print(*output, sep="\n", file=tempout)
        os.chmod(args.script, 0o755)
        log_message(f"created {args.script}")
    else:
        print(*output, sep="\n")


@cyc_app.command(group=cyc_group_refs)
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
        log_message(f"Mismatch in DNA and RNA files for {species} aka {species_Latin}:\nDNA: {dna}\nRNA{rna}", fatal=True)

    dna_build = os.path.basename(dna).split(".")[1]
    dna_build_from_rna_file, rna_build = os.path.basename(rna).split(".")[1:3]

    if dna_build == dna_build_from_rna_file:
        log_message(f"DNA and RNA files are both for genome build {dna_build}")
    else:
        log_message(f"Genome builds mismatched {dna_build} vs. {dna_build_from_rna_file}:\nDNA: {dna}\nRNA{rna}", fatal=True)
    
    ens = []
    for f in [dna, rna]:
        if m := re.search(r'release-(\d+)\/', f):
            ens.append(m.groups()[0])
        else:
            log_message(f"Release not found in {f}", fatal=True)
    if len(set(ens)) != 1:
        log_message(f"Ensembl releases mismatched {ens[0]} vs. {ens[1]}:\nDNA: {dna}\nRNA{rna}", fatal=True)
    ens = ens[0]
    return dna_build, ens
    




@cyc_app.command(group=cyc_group_align)
def bam_species(*, bamfiles: list[FileName], outputfile: FileName | None = None):
    if outputfile:
        exit_if_files_exist(outputfile)
    verify_that_paths_exist(bamfiles)
    species = get_species_for_bamfiles(bamfiles=bamfiles)
    with open(outputfile, "w") if outputfile else sys.stdout as tempout:
        for bam in bamfiles:
            print(f"{bam}\t{species[bam]}", file=tempout)

#check-if-ref-files-match
@cyc_app.command(group=cyc_group_align)
def bash_header(*, flags: str = "set -e"):
    return dedent(f"""
        #!/bin/bash
        {flags}
        log() {{ echo "$*" >&2 ; }}
        die() {{ echo "$*" >&2; exit 1 ; }}
    """).lstrip()


@dataclass
class get_genecoords:
    genes: list[str]
    "Gene(s) in a file or as a space-separated list"

    species: Literal[constants.known_species]
    "On the origin of reads."

    extend: int = 0
    "Extension of gene boundaries, in nucleotides. Maximum not checked against chromosome length."

    outputfile: File_or_Dir = "gene_coords_{species}_ext_{extend}.bed"
    "Output file for gene coordinates. Use double quotes for stdout."

    overwrite: bool = False
    "Overwrite the output file if it's already present."

    script: str | None = None
    'Output command to this Bash script instead of executing.'

@cyc_app.command(group=cyc_group_refs)
def get_genecoords(args):
    """Get coordinates for a list of genes, in BED format, or output the script to do this."""


    #if args.script and args.outputfile:
    #    log_message("Specify either --script or --outputfile but not both.", fatal=True)
    outputfile = args.outputfile
    outputfile = outputfile.replace("{species}", args.species)
    outputfile = outputfile.replace("{extend}", str(args.extend))
    if not args.overwrite:
        if outputfile:
            exit_if_files_exist(outputfile)
        if args.script:
            exit_if_files_exist(args.script)
    species = args.species
    coords_file = os.path.expandvars(constants.coords_source_file[species])
    verify_that_paths_exist(coords_file)
    genes_to_get = return_list_from_file_or_list(args.genes)
    if args.script:
        genes = " ".join(genes_to_get)
        output = f"""
            #!/bin/bash
            set -e
            ext={args.extend}
            genes={genes}
            species={species}
            {__file__} get_genecoords -g $genes -s $species -x $ext -o {outputfile}
        """
        output = dedent(output)
        with open(script, "w") as tempout:
            print(output, file=tempout)
        os.chmod(script, 0o755)
        log_message(f"created {args.script}")
    else:
        lowercase = set(x.lower() for x in genes_to_get)

        log_message(f"reading {species} gene coordinates from {coords_file}")
        data = pd.read_csv(coords_file, sep="\t", header=0)
        if " ".join(data.columns[0:4]) != "chr start stop gene_name":
            log_message(
                f"Unexpected file contents in  {coords_file}", fatal=True
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
        data["start"] = data["start"].astype(int)
        data["stop"] = data["stop"].astype(int)
        if args.extend:
            data["start"] -= args.extend
            data["start"] = [
                max(1, x) for x in data["start"]
            ]  # data.apply(lambda x: max(x["start"], 1), axis=1)
            data["stop"] += args.extend
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

@cyc_app.command(group=cyc_group_align)
def rnaspades(args):

    # construct commands to:
    # 1. extract unaligned reads from starting BAM (i.e., reads not mapped as proper pairs)
    # 2. extract all paired reads from extended gene-specific BAM, discarding singletons (from ends)
    # 3. run rnaspades
    # Only one fastq file is generated for each sample, but it matters for rnaspades whether the reads are single- or paired-end.
    outputscript = args.script
    if outputscript:
        exit_if_files_exist(outputscript)
    test_executables(constants.rnaspades.exe)

    genebamfiles = args.genebamfiles or []
    genebamfiles = [x for x in genebamfiles if is_a_gene_bam(x)]
    bamfiles = args.bamfiles or []
    bamfiles = [x for x in bamfiles if not is_a_gene_bam(x)]
    regions_file = args.regions or ""

    if bamfiles and genebamfiles:
        genebam = {x: gene_bam(x) for x in bamfiles}
        if set(genebam.values()) == set(genebamfiles):
            log_message("Standard and gene BAM files match")
        else:
            log_message(
                "Problem - standard and gene BAM files specified but mismatched", fatal=True
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
                    "Specify gene bams or region file to create them", fatal=True
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
    for d in outputdir_by_sample.values(): check_dir_write_acycess(d)

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
    cpus = args.cpus
    mem = args.mem

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
        if args.zipbams:
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
        execute =args.exec,
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
    bamfiles = [x for x in args.bamfiles if is_a_gene_bam(x) == False]
    verify_that_paths_exist(bamfiles + [f"{bam}.bai" for bam in bamfiles])
    outputscript = args.script
    if outputscript:
        exit_if_files_exist(outputscript)
    regions_file = args.regions
    verify_that_paths_exist(regions_file)
    outputbam = {x: gene_bam(x) for x in bamfiles}
    exit_if_files_exist(list(outputbam.values()))
    output = []
    if args.exec:
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

    if args.zipbams:
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
            execute =args.exec,
        )
    # else:
    #    log_message("Error - nothing to do", fatal=True)


def junctions(args):
    pass
    from tempfile import mkstemp

    # needs coords and either gene bams or non-gene BAMs, but not both
    test_executables("regtools")
    outputfileext = args.ext
    bamfiles = args.bamfiles
    verify_that_paths_exist(bamfiles)
    outputfile = {bamfile: f"{bamfile}.{outputfileext}" for bamfile in bamfiles}
    # outputfileprefix = {x: re.sub('\.bam$', f".realigned", x) for x in genebamfiles}
    exit_if_files_exist(outputfile.values())
    mincount = args.mincount
    if len(set(outputfile.values())) != len(set(outputfile.keys())):
        log_message("Redundant output file(s)", fatal=True)
    regions = args.regions
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
    genebamfiles = args.genebamfiles or []
    genebamfiles = [x for x in genebamfiles if is_a_gene_bam(x)]
    bamfiles = args.bamfiles or []
    bamfiles = [x for x in bamfiles if not is_a_gene_bam(x)]

    outputscript = args.script or ""
    if outputscript:
        exit_if_files_exist(outputscript)
    regions_file = args.regions
    verify_that_paths_exist(regions_file)
    cpus = args.cpus - 1
    if bamfiles:
        genebam_for_bam = {x: gene_bam(x) for x in bamfiles}
        if genebamfiles:
            if set(genebam_for_bam.values()) == set(genebamfiles):
                log_message("Standard BAMs ignored since gene bams exist")
                bamfiles = []
            else:
                log_message(
                    "Problem - bams & gene bams specified but mismatched", fatal=True
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
        log_message("No input BAM or gene BAM file(s)", fatal=True)
    verify_that_paths_exist(genebamfiles)
    outputbam = {x: realigned_gene_bam(x) for x in genebamfiles}
    exit_if_files_exist(list(outputbam.values()))
    species = args.species or get_single_species_for_bamfiles(bamfiles=genebamfiles)
    read_type = args.reads or get_single_read_type_for_bamfiles(bamfiles=genebamfiles)
    command = f"abra2 --nosort --threads {cpus} --undup"
    if read_type == constants.read_type.single:
        command += " --single"
    elif read_type != constants.read_type.paired:
        log_message(f"Invalid read type {read_type}", fatal=True)
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
        log_message(f"ref files not found:\n{errors}", fatal=True)

    outputfileprefix = {x: re.sub(r"\.bam$", f".realigned", x) for x in genebamfiles}
    temp_region_file_by_chr = {chr: f"{regions_file}.temp.{chr}" for chr in chrs}
    if len(set(temp_region_file_by_chr.values())) != len(
        set(temp_region_file_by_chr.keys())
    ):
        log_message(f"length mismatch {str(temp_region_file_by_chr)}", fatal=True)

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
        if args.zipbams:
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
        execute =args.exec,
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
            log_message(f"Fewer than two columns found in {file}", fatal=True)
        return (
            temp[temp.columns[0:2]].set_index(temp.columns[0]).T.to_dict("records")[0]
        )
    # log_message(f"Obtaining {arglabel} from command line")
    # expecting from:to pairs
    rename = {}
    for pair in argvalues:
        temp = pair.split(":")
        if len(temp) != 2:
            log_message(f"Malformed from:to pair {pair}", fatal=True)
        rename[temp[0]] = temp[1]
    return rename


def find_gtf(*, species : str):
    gtf = constants.gtf.get(species, "")
    if not gtf:
        log_message(f"Undefined default gtf for {species}", fatal=True)
    if os.path.exists(gtf):
        return gtf
    if os.path.exists(f"{gtf}.gz"):
        log_message(f"""
            zcat {gtf}.gz > {gtf}
               or
            gunzip {gtf}
        """, fatal=True)
    log_message(f"{gtf} not found for {species}", fatal=True)


def featureCounts(args):
    # expects all BAM files to be from the same species
    outputscript = args.script
    if outputscript:
        exit_if_files_exist(outputscript)
    bamfiles = args.bamfiles
    verify_that_paths_exist(bamfiles)

    sort_by_read_ID = args.sort_BY_READ_ID == "yes"
    if sort_by_read_ID:
        unsortedbam = {bam: re.sub(r"\.bam$", ".unsorted.bam", bam) for bam in bamfiles}
        # already_sorted = [bam for bam in bamfiles if os.path.exists(unsortedbam[bam])]
        outputfile = {
            bam: re.sub(r"\.bam$", ".unsorted.counts.txt", bam) for bam in bamfiles
        }
    else:
        outputfile = {bam: re.sub(r"\.bam$", ".counts.txt", bam) for bam in bamfiles}
    exit_if_files_exist(outputfile.values())
    species = args.species or get_single_species_for_bamfiles(bamfiles=bamfiles)

    # constants.sortcpus = floor(args.cpus / 2)

    
    # log_message(f"constants.sortcpus: {constants.sortcpus}\nconstants.sortmem: {constants.sortmem}")
    if (args.sort_CPUS + 1) * args.sortmem > system_memory():
        log_message(
            f"{args.sort_CPUS +1} CPUs x {args.sortmem}M > {system_memory()}", fatal=True
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
                samtools view -F 0x4 --bam {bam} 2> {out}.unsort.1.err | samtools sort -n -@ {args.sort_CPUS} -m {args.sortmem}M -l 9 -o {tempbam} 2> {out}.unsort.2.err
                featureCounts {constants.featurecounts.options} -a {gtf} -o {out} -T {args.cpus} {tempbam} 2> {out}.log
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
        execute =args.exec,
    )


def unsort(args):
    outputscript = args.script
    if outputscript:
        exit_if_files_exist(outputscript)
    bamfiles = args.bamfiles
    verify_that_paths_exist(bamfiles)
    outputfile = {bam: re.sub(r"\.bam$", ".unsorted.bam", bam) for bam in bamfiles}
    exit_if_files_exist(outputfile.values())

    sortcpus = args.sort_CPUS
    sortmem = args.sortmem
    log_message(
        f"sortcpus: {sortcpus}\nsortmem: {sortmem}"
    )
    if (sortcpus + 1) * sortmem > system_memory():
        log_message(
            f"{sortcpus +1} CPUs x {sortmem}M > {system_memory()}",
            fatal=True,
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
        execute =args.exec,
    )


def unique_gene_metadata(x):
    return ";".join(list(set(str(x).split(";"))))

def gene_metadata(args):
    log_message("output gene metadata")
    # inputfile = args.inputfiles[0]
    verify_that_paths_exist(args.inputfile)
    outputfile = args.outputfile
    if outputfile:
        if (b := os.path.basename(outputfile)) != outputfile:
            log_message(f"outputfile {outputfile} != basename {b}", fatal=True)
        outputdir = args.outputdir.rstrip("/")
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
        log_message(f"Insufficient columns in {inputfile}", fatal=True)
    if len(data.columns) > 3:
        log_message(f"Too many columns found in {inputfile}", fatal=True)
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
        log_message("Invalid number of columns", fatal=True)
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
        log_message("Errors found in fastq file manifest", fatal=True)

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


@cyc_app.command(group=cyc_group_refs)
def star_make_idx(args):
    # Output command to make an index for a given species and read length.
    if args.script:
        args.script = args.script.replace("{species}", args.species)
        args.script = args.script.replace("{readlength}", str(args.readlength))
        if not args.overwrite:
            exit_if_files_exist(args.script)
    index = args.index
    index = index.replace("{species}", args.species)
    index = index.replace("{readlength}", str(args.readlength))

    dna = constants.ref.dna[args.species]
    rna = constants.ref.rna[args.species]
    print(f"dna={dna}, rna={rna}, species={species}")
    dna_build, rna_build = check_if_ref_files_match(dna=dna, rna=rna, species=args.species)
    #"human": "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz",
    #"mouse": "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
    dna = os.path.basename(dna).replace(".gz","")
    rna = os.path.basename(rna).replace(".gz","")
    index = index.replace("{dna_build}", dna_build)
    index = index.replace("{rna_build}", rna_build)


    cpus = args.cpus or "$(grep -c ^processor /proc/cpuinfo)"
    if args.mem == 0:
        mem=""
    else:
        mem = f"--limitGenomeGenerateRAM {args.mem}"
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
            --genomeFastaFiles $dna --sjdbGTFfile $rna --sjdbOverhang {args.readlength} \\
            --runThreadN $cpus {mem} > $index.log 2> $index.err
    '''))
    if args.script:
        with open(args.script, "w") as tempout:
            print("\n".join(output), file=tempout)
        os.chmod(args.script, 0o755)
        log_message(f"created {args.script}")
    else:
        print("\n".join(output))

    _

@cyc_app.command(group=cyc_group_align)
def star(args):
    # output commands to fetch fastqs and run STAR and featureCounts
    test_executables(["samtools", "STAR"], exit_on_error=True)
    verify_that_paths_exist(
        [constants.star.base_dir, constants.bin_dir, constants.ref_dir],
        exit_on_error=True,
    )
    cpus = args.cpus
    sortcpus = args.sort_CPUS
    sortmem = args.sortmem
    # prefetch = args.prefetch
    verify_that_paths_exist(args.inputfile)
    if args.abra2COORDS:
        verify_that_paths_exist(args.abra2COORDS)
        test_executables("abra2", exit_on_error=True)
    """
    if args.rnaspadescoords:
        verify_that_paths_exist(args.rnaspadescoords)
        test_executables(rnaspades.exe)
    """
    outputdir = args.outputdir
    counts_dir = os.path.join(outputdir, "counts")
    check_dir_write_access([outputdir, counts_dir])
    if glob.glob(f"{counts_dir}/*"):
        log_message(f"\nCaution: files found in {counts_dir}\n", fatal=True)
    bamdestdir = args.bamdestdir or ""
    if bamdestdir:
        check_dir_write_access(bamdestdir)
    if args.abra2COORDS:
        abra2dir = os.path.join(outputdir, "abra2")
        check_dir_write_access(abra2dir)
    """
    if args.rnaspadescoords:
        rnaspadesdir = os.path.join(outputdir, "rnaspades")
        check_dir_write_access(rnaspadesdir)
        # rnaspadesdir = os.path.join(outputdir, "rnaspades")
        # temp = os.path.basename(os.path.realpath(outputdir))
        # rnaspadesdir = check_dir_write_access(dir = os.path.join(rnaspadesdir, temp))
    """
    outputscript = args.script  # os.path.join(outputdir, )
    overwrite = args.overwrite
    if outputscript and not overwrite:
        exit_if_files_exist(outputscript)
    indexpath, species = find_star_index(args.index)
    # find gtf
    gtf = find_gtf(species=species)
    fastq_source = args.fastqsource.upper()
    samples, fastq_fetch_params, star_fastq_files, read_type = parse_fastq_manifest(
        inputfile=args.inputfiles,
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
        function_name, function_code = bash_aspera(speed=args.transfer_SPEED)
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
        cpus=args.cpus, index=indexpath, readFilesCommand=readFilesCommand
    )
    bash_functions["star"] = function_name
    output.append(function_code)

    function_name, function_code = bash_samtools_sort_by_position(
        sortcpus=sortcpus, sortmem=sortmem
    )
    bash_functions["sort BAM by read ID"] = function_name
    output.append(function_code)

    if args.counts:
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
        if args.add_OPTS:
            other_opts.append(args.add_OPTS)
        if args.addreadgroup:
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
        if args.counts:
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

    if bamdestdir and args.abra2COORDS: # or args.rnaspadescoords):
        output.append(dedent(f"""
            log waiting for mvjobid=$mvjobid
            wait $mvjobid
            log $mvjobid done
            # temp = os.path.realpath(bamdestdir)
        """))

    if args.abra2COORDS:
        srcdir = bamdestdir or os.path.realpath(outputdir)
        output.append(f"cp -s {srcdir}/*.bam {srcdir}/*.bai {abra2dir}/")
        output.append(f"mv {args.abra2COORDS} {abra2dir}/")
        output.append(f"cd {abra2dir}")
        output.append(
            f"{__file__} abra2 -b *.bam -r {os.path.basename(args.abra2COORDS)} -z --exec"
        )
        output.append(f"cd {os.path.realpath(outputdir)}")
    """
    if args.rnaspadescoords:
        srcdir = bamdestdir or os.path.realpath(outputdir)
        output.append(f"cp -s {srcdir}/*.bam {srcdir}/*.bai {rnaspadesdir}")
        output.append(f"mv {args.rnaspadescoords} {rnaspadesdir}")
        output.append(f"cd {rnaspadesdir}")
        output.append(
            f"{__file__} rnaspades -b *.bam -r {os.path.basename(args.rnaspadescoords)} -z --exec"
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
    files = args.inputfiles
    verify_that_paths_exist(files)
    outputfile = args.outputfile
    if outputfile:
        outputdir = args.outputdir.rstrip("/")
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
        log_message("\nCaution: non-unique sample IDs\n", fatal=True)

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
            f"{outputfile} exists - delete it or allow overwrite", fatal=True
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
    verify_that_paths_exist(args.inputfile)
    outputfile = args.outputfile
    # exit_if_files_exist(outputfile)
    # header = None if args.noheader else 0
    data = pd.read_csv(args.inputfile, sep="\t", header=0)
    sum_counts(
        data=data,
        sum_column="total_counts",
        outputfile=args.outputfile,
        overwrite=args.overwrite,
    )
    log_message(f"Total counts written to {outputfile or 'stdout'}")


def add_gene_name(args):
    log_message("substitute gene name in expression counts if unique")
    verify_that_paths_exist(args.inputfile)
    outputfile = args.outputfile
    if outputfile:
        exit_if_files_exist(outputfile)
    species = args.species
    ref_file = f"${HOME}/star/ref/ens.to.gene.{species}.if.unique"
    gene_names = pd.read_csv(
        ref_file, sep="\t", header=None, names=["Geneid", "dummy"], dtype=object
    )
    data = pd.read_csv(args.inputfile, sep="\t", header=0, dtype=object)
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
    files = args.inputfiles
    verify_that_paths_exist(files)
    outputfile = args.outputfile
    outputdir = args.outputdir.rstrip("/")
    totals = args.totals
    if outputdir:
        #if outputfile and os.sep in outputfile:
        if outputfile and os.path.basename(args.outputfile) != args.outputfile:
            log_message(
                "outputfile cannot include dir if outputdir is specified", fatal=True
            )
        if totals and os.path.basename(totals) != totals:
            log_message(
                "output file for totals cannot include dir if outputdir is specified",
                fatal=True,
            )
        totals = args.totals
        outputfile = os.path.join(outputdir, outputfile)
        totals = os.path.join(outputdir, totals)
        exit_if_files_exist([outputfile, totals])
    output = get_one_file_of_counts(file=files[0])
    for file in files[1:]:
        data = get_one_file_of_counts(file=file)
        if common_samples := set(data.columns[1:]) & set(output.columns[1:]):
            common_samples = "\n".join(common_samples)
            log_message(f"Common sample IDs:\n{common_samples}", fatal=True)
        output = output.merge(data, how="inner", on="Geneid")
    rename_samples = args.rename_SAMPLES
    if rename_samples:
        rename_samples = pd.read_csv(
            args.rename_SAMPLES,
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
    verify_that_paths_exist(args.inputfile)
    if args.outputfile and not args.overwrite:
        exit_if_files_exist(args.outputfile)

    temp = pd.read_csv(args.inputfile, sep="\t", header=0, comment="#", nrows=1)
    if args.column:
        verify_columns_present(
            data=temp, columns=args.column, source=args.inputfile
        )
        if not args.newcolumn:
            newcolumn = f"rank_by_{args.column}"
    else:
        if args.sort:
            log_message("There is no sense in sorting by the rank", fatal=True)
        newcolumn = args.newcolumn or "order"
    verify_columns_absent(
            data=temp,  columns=newcolumn, source=args.inputfile
        )
    data = pd.read_csv(args.inputfile, sep="\t", header=0, comment="#")

    if args.column:
        if args.sort:
            data.sort_values(by=args.column, ascending = args.order == "ascending", inplace=True)
        data[newcolumn] = data[args.column].rank(method="min").astype(int)
    else:
        data[newcolumn] = list(range(1, data.shape[0]+1))
    data.to_csv(args.outputfile if args.outputfile else sys.stdout, sep="\t", index=False)


def transformcounts(args):
    log_message("transform counts")
    verify_that_paths_exist(args.inputfile)
    method = args.method.lower()
    if "log2" in method and args.log2_OFFSET <= 0:
        log_message(
            "log2_offset must be positive when applying log", fatal=True
        )

    outputfile = args.outputfile
    if not outputfile:
        outputfile = f"{method.lower()}.txt"
    overwrite = args.overwrite
    if outputfile and overwrite is False:
        exit_if_files_exist(outputfile)

    # get metadata
    metadata = args.metadata
    if metadata:
        metadata = os.path.join(os.path.dirname(args.inputfile), metadata)
        verify_that_paths_exist(metadata)
        use_genenames = args.gene_NAMES if metadata else 0
        metadata_columns_to_get = ["Geneid", "Length", "gene_biotype"]
        if use_genenames:
            metadata_columns_to_get.append("gene_name")
        temp = pd.read_csv(metadata, sep="\t", header=0, comment="#", nrows=1)
        verify_columns_present(
            data=temp, columns=metadata_columns_to_get, source=metadata
        )
        metadata = pd.read_csv(
            metadata, sep="\t", header=0, usecols=metadata_columns_to_get, comment="#"
        )
        gene_types = args.gene_TYPES
        if gene_types and gene_types != "all":
            mask = metadata["gene_biotype"].isin(gene_types)
            metadata = metadata[mask].copy()
    data = pd.read_csv(args.inputfile, sep="\t", header=0)
    #  sum counts in case gene IDs are repeaed
    data = data.groupby(data.columns[0], as_index=False).sum()
    samples = [
        x for x in data.columns[1:] if not x in constants.featurecounts_constant_columns
    ]

    # ["CPM", "CPM-UQ", "CPM-UQ-log2", "RPKM", "RPKM-UQ", "RPKM-UQ-log2"]
    temp_output = f"{args.totals}.txt"
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

        temp_output = f"{args.totals}_filtered.txt"
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
                log_message(f"Unknown method {method}", fatal=True)

        temp_output = f"{args.totals}_sum_{method}.txt"
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
        pct = np.nanpercentile(tempdata, args.rescale_PERCENTILE, axis=0)
        data[samples] = data[samples] * args.rescale_COMMON_VALUE / pct

    if "log2" in method:
        data[samples] = np.log2(data[samples].applymap(lambda x: x + args.log2_OFFSET))

    data.round(2).to_csv(outputfile or sys.stdout, sep="\t", index=False)

    temp = outputfile or "stdout"
    log_message(f"{method} written to {temp}")


def counts_postproc(args):
    output = []
    if args.script:
        exit_if_files_exist(args.script)
    dir = args.dir.rstrip("/")
    species = args.species
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
    cmd = f"{args.script} &>> {args.script}.log" if args.script else ""
    handle_commands_and_output(
        commands=output,
        outputfile=args.script,
        single_line_command=cmd,
        execute =args.exec,
    )

@cyc_app.command(group=cyc_group_align)
def star_zip_files(args):
    dir = args.dir.rstrip("/")
    outputscript = args.script
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
        execute =args.exec,
    )

cyc_group_align = cyclopts.Group.create_ordered("")
def star_clear_mem(args):
    outputscript = args.script
    if outputscript:
        exit_if_files_exist(outputscript)
    verify_that_paths_exist(constants.star.dummy_fastq)
    output = [bash_header()]
    output.append(dedent(f"""
        STAR --genomeDir {args.index} --readFilesIn {constants.star.dummy_fastq} --runThreadN 4 --outFile_or_DirPrefix ${constants.star.dummy_fastq}. --genomeLoad Remove --runMode alignReads"
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
    check_for_file_xor_stdin(args.inputfile)
    outputfile = args.outputfile
    if outputfile:
        exit_if_files_exist(outputfile)
    # seq = ""
    minlength = args.min_LENGTH
    with open(args.inputfile, "r") if args.inputfile else sys.stdin as tempin:
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

@dataclass
class fasta_to_text:
    inputfile: File_or_Dir | None = None
    "Input file."

    minlength: int = 0
    "Minimum length"

    outputfile: File_or_Dir | None = None
    "Output file. Stdout if omitted."
# join

@cyc_app.command()
def fasta_to_text(args):
    help_text ="Convert fasta file to tab-delimited text."
    # inputs fasta file, outputs text file with header row
    # presumably called from command line
    # output to a file
    if args.inputfile:
        verify_that_paths_exist(args.inputfile)
    else:
        detect_stdin(exit_on_error=True)
    outputfile = args.outputfile
    if outputfile:
        exit_if_files_exist(outputfile)
    minlength = args.min_LENGTH
    seq = ""
    if minlength:
        log_message(f"Caution: min length = {minlength}")
    with open(args.inputfile, "r") if args.inputfile else sys.stdin as tempin:
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




"""
def join(args):
    inputfiles = args.inputfiles
    if len(inputfiles) != 2:
        log_message("Specify two input files", fatal=True)
    verify_that_paths_exist(inputfiles)
    outputfile = args.outputfile
    if outputfile: exit_if_files_exist(outputfile)
    method = args.method
    columns = args.columns
    if len(columns) == 2:
        for i, cols in enumerate(columns):
            columns[i] = cols.split(",")
    elif len(columns) == 1:
        columns.append(columns[0].split(","))
        columns[0] = columns[0].split(",")
        #columns += columns[0]
    else:
        log_message("Invalid columns", fatal=True)
    data = []
    for i in [0, 1]:
        data.append(pd.read_csv(inputfiles[i], sep="\t", header = 0, dtype=object))
        #data[i].to_csv(f"out.{i}", sep="\t", index=False)
        #print(f"{i}:\n" + "\n".join(data[i].columns))
        verify_columns_present(data = data[i], columns=columns[i], source = inputfiles[i])
    output = pd.merge(data[0], data[1], left_on = columns[0], right_on=columns[1], how = method)
    if output.empty:
        log_message("No common column values", fatal=True)
    if args.dedup:
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
    # outputfileext = args.ext
    verify_that_paths_exist(args.inputfiles)
    outputscript = args.script
    exit_if_files_exist(outputscript)
    species = args.species
    junctions = constants.minimap_junctions[species]
    genome = constants.minimap_genome[species]
    if genome == "":
        log_message(f"No minimap genome index for {species}", fatal=True)
    if junctions == "":
        log_message(f"No minimap junctions for {species}", fatal=True)
    genome = os.path.join(constants.ref_dir, genome)
    junctions = os.path.join(constants.ref_dir, junctions)
    verify_that_paths_exist([genome, junctions])
    outputformat = args.outputformat
    # outputfileprefix = {x: minimap_output_file(inputfile=x, ext=outputfileext) for x in inputfiles}
    outputfileprefix = {x: minimap_output_file(inputfile=x) for x in args.inputfiles}
    exit_if_files_exist([f"{x}.{outputformat}" for x in outputfileprefix.values()])
    output = []
    output.append(bash_header())
    function_name, function_code = bash_minimap(
        cpus=args.cpus, bases_loaded=args.mem, genome=genome, junctions=junctions
    )
    output.append(function_code)
    for inputfile in args.inputfiles:
        output.append(f"# {inputfile}")
        output.append(
            f"{function_name} {inputfile} {outputfileprefix[inputfile]} {outputformat}"
        )
    cmd = f"nohup {outputscript} &>> {outputscript}.log &"
    handle_commands_and_output(
        commands=output,
        outputfile=outputscript,
        single_line_command=cmd,
        execute =args.exec,
    )

    # bam2Bed12.py -i transcripts.fasta.v4.bam > transcripts.fasta.v4.bed2.txt
'''

def ens_to_gene(args):
    gtf_chr =  0
    gtf_feature_type =  2
    gtf_features =  -1
    if args.species:
        if args.inputfile or args.outputfile:
            log_message("Specify species, or specify input and output files.", fatal=True)
        inputfile = os.path.expandvars(constants.gtf[args.species])
        if not os.path.exists(inputfile):
            if os.path.exists(inputfile + ".gz"):
                cmd = f"gunzip --keep {inputfile}.gz"
                log_message(f"Executing {cmd}")
                execute_command(cmd)
        outputfile = constants.featurecounts.unique_IDs[args.species]
        outputfile = os.path.expandvars(outputfile)
    else:
        if args.inputfile is None or args.outputfile is None:
            log_message("Specify both input file and output file, or specify species.", fatal=True)
        inputfile = args.inputfile
        outputfile = args.outputfile
    verify_that_paths_exist(inputfile)
    if not args.overwrite:
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
    if args.species:
        if args.inputfile or args.outputfile:
            log_message("Specify species, or specify input and output files.", fatal=True)
        inputfile = os.path.expandvars(constants.gtf[args.species])
        if not os.path.exists(inputfile):
            if os.path.exists(inputfile + ".gz"):
                cmd = f"gunzip --keep {inputfile}.gz"
                log_message(f"Executing {cmd}")
                execute_command(cmd)
        outputfile = constants.coords_source_file[args.species]
        outputfile = os.path.expandvars(outputfile)
    else:
        if args.inputfile is None or args.outputfile is None:
            log_message("Specify both input file and output file, or specify species.", fatal=True)
        inputfile = args.inputfile
        outputfile = args.outputfile
    verify_that_paths_exist(inputfile)
    if not args.overwrite:
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
    if args.genes:
        genes_before = data.shape[0]
        biotypes = data["gene_biotype"].unique().tolist()
        if common := set(args.genes) & set(biotypes):
            if missing := set(args.genes) - set(biotypes):
                log_message("Biotypes not found in GTF:\n\t" + "\n\t".join(list(missing)))
            mask = data["gene_biotype"].isin(common)
            data = data[mask].copy()
        else:
            log_message("None of the biotypes are not found in the GTF.", fatal=True)
        log_message(f"{data.shape[0]} of {genes_before} genes retained")
    if args.sort:
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
        log_message(f"Unknown species {species}", fatal=True)
    for dir in constants.index_stores:
        indexpath = os.path.join(dir, index)
        if os.path.exists(indexpath):
            if dir != constants.star.base_dir:
                log_message(
                    f"nice cp -r {indexpath} {constants.star.base_dir} &"
                )  # , fatal=True)
            return indexpath, species
    temp = "\n".join(constants.index_stores)
    log_message(f"{index} not found in:\n{temp}", fatal=True)

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





if __name__ == "__main__":
    constants.star_indexes = find_star_indexes()
    cyc_app()


