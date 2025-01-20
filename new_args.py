
# work in progress - not implemented
# This is how I would implement command-line arguments today.

from cyclopts import App, Group

app = App()

metadata = Group.create_ordered("metadata")
ref = Group.create_ordered("ref")
align = Group.create_ordered("align")
counts = Group.create_ordered("counts")
assembly = Group.create_ordered("assembly")
bamfiles = Group.create_ordered("bamfiles")
utils = Group.create_ordered("utils")


@app.command(group=metadata)
help_text ="Find metadata and files available from GEO, SRA and ENA."

@dataclass
class metadata:


    studydir: FileName
    "Folder where metadata files will be stored, with each source in its own subfolder."


    ID: str
    "GEO, SRA or ENA identifier. Ex. PRJNA743892 GSE162198 SRP294329"


    sources: Union[Literal[*constants.metadata_sources]] = # constants.metadata_sources # = list[str]
    "Source(s) of metadata to query"


    assay: str | list[str]  =  ["RNA-Seq"]
    "Data type of interest."


    species: Literal[constants.known_species]
    "Species of interest."


# metadata

@app.command(group=metadata)
@dataclass
class pysradb:
help_text ="Get metadata with pysradb"


    IDs: list[str]
    "Identifier(s)."


    keeptempfiles: bool = False
    "Keep temp output files."


    outputfile: FileName | None = None
    "Output file for a single ID. Default is {ID}.{ext}"


    ext: str = "_pysradb.txt"
    "Output file suffix"


    overwrite: bool = False
    "Overwrite files if present."

@app.command(group=metadata)
@dataclass
class geo:
help_text ="Get metadata from GEO"


    ID: list[str]
    "Identifier(s)."


    ext: str = "_geo.txt"
    "Output file suffix"

# dedup_cols
@app.command(group=metadata)
@dataclass
class dedup_cols:
help_text ="Remove duplicate columns"

    inputfile: FileName, | None = None
    "Input file."


    outputfile: FileName, | None = None
    "Output file. Stdout if omitted."
# enafastqs


@app.command(group=metadata)
@dataclass
class enafastqs:
help_text ="Get a ist of ENA fastqs for ID(s) like PRJNA627881"


    IDs: list[str], required=True
    "Identifier(s)."


    keeptempfiles: bool = False
    "Keep temp output files."


    ext: str = "ena.txt"
    "Output file extension"

# genecoords
@app.command(group=metadata)
@dataclass
class get_genecoords:
help_text ="Get coordinates for a list of genes, in BED format, or output the script to do this."


    genes: list[str]
    "Gene(s) in a file or as a space-separated list"


    species: Literal[constants.known_species]
    "On the origin of reads."


    extend: int = 0
    "Extension of gene boundaries, in nucleotides. Maximum not checked against chromosome length."


    outputfile: FileName = "gene_coords_{species}_ext_{extend}.bed"
    "Output file for gene coordinates. Use double quotes for stdout."


    overwrite: bool = False
    "Overwrite the output file if it's already present."

    script: str | None = None
    'Output command to this Bash script instead of executing.'


# gene_metadata
@app.command(group=count)
@dataclass
class gene_metadata:
help_text ="Extract metadata columns from featureCounts output files."

    inputfile: FileName
    "Input file."


    outputdir: str = "."
    "Output directory."


    outputfile: FileName = constants.default_gene_metadata
    "Output file for gene metadata."

# pre-STAR
@app.command(group=ref)
@dataclass
class genomeref:
help_text ="Generate script to fetch genome (fasta) and transcriptome (gtf) files."


    species: Literal[constants.known_species]
    "On the origin of reads."


    outputdir: FileName = constants.ref_dir
    "Output directory."


    overwrite: bool = False
    "Overwrite files if present."

    script: str | None = f"genomeref_{{species}}.sh"
    'Output commands to a Bash script. Set to "" for stdout.'

@app.command(group=ref)
@dataclass
class gtf_to_coords:
help_text ="Extract gene coordinates from GTF."

    inputfile: FileName, | None = None
    "Input file. Required unless --species is specified."


    outputfile: FileName, | None = None
    "Output file.  Required unless --species is specified."


    species: str | None = None, Literal[constants.known_species]
    "On the origin of reads. Assumes default input & output files."


    overwrite: bool = False
    "Overwrite files if present."


    sort : bool = True
    "Sort by chr, start & stop. Use --no_sort to skip."


    genes: list[str] = ["protein_coding"]
    "Gene type(s) to keep, space-separated. Use \"\" for all."

@app.command(group=ref)
@dataclass
class ens_to_gene:
help_text ="Extract unique gene-name/gene-ID pairs from GTF."

    inputfile: FileName, | None = None
    "Input file. Required unless --species is specified."


    outputfile: FileName, | None = None
    "Output file.  Required unless --species is specified."


    species: str | None = None, Literal[constants.known_species]
    "On the origin of reads. Assumes default input & output files."


    overwrite: bool = False
    "Overwrite files if present."

# STAR alignments etc.
@app.command(group=align)
@dataclass
class star:
help_text ="Output commands to fetch fastqs then run STAR & featureCounts. Input=label tab fastq1 (tab fastq2)"


    index: str required=True, Literal[constants.star_indexes]
    "Genome index. Path is ignored."

    inputfile: FileName
    "Input file."


    abra2: FileName | None = None
    "Coords file to use for gene bams then abra2."


    additional: str = ""
    "Additional options, in double quotes, added verbatim."


    addreadgroup: bool = False
    "Use sample ID as read group."


    bamdest: FileName | None = None
    "Destination for BAM files after completion, e.g. on mounted drive"


    counts: bool = True
    "Run featurecounts after alignments. Use --no_counts to skip."


    cpus: int = constants.cpus
    "CPUs."


    exec: bool = False
    "Execute command(s) instead of output."


    fastqsource: Literal["ENA", "SRA"] = "ENA"
    "Source of fastq files."


    outputdir: FileName = "."
    "Output directory."


    overwrite: bool = False
    "Overwrite files if present."

    readtype: Literal[constants.read_types] | None = None
    "Read pairing"


    rnaspades: FileName | None = None
    "Coords file to use for gene bams then rnaspades."


    sortcpus: int = int(constants.cpus / 2)
    "CPUs to use for sorting"

    sortmem: int  = 2000
    "Memory to use for sorting"


    transferspeed: int = 500
    "Aspera transfer speed."


    script: FileName = default=f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'

@app.command(group=align)
@dataclass
class star_clear_mem:
help_text ="Output command to remove star genome from mem"

    index: str required=True
    "Genome index. Path is ignored."


    script: FileName = "star_clear_mem.sh"
    'Output commands to bash script. Set to "" for stdout.'

@app.command(group=align)
@dataclass
class star_make_idx:
help_text ="Output command to make an index for a given species and read length."


    species: Literal[constants.known_species]
    "On the origin of reads."


    readlength: int
    "Read length minus 1, probably 49, 75, 100 or  150. See Star Manual"


    index: str = os.path.join(constants.star.base_dir, "{species}.{dna_build}_{rna_build}.{readlength}")
    "Name of index, i.e., base output directory."


    cpus: int = constants.cpus
    "CPUs. Set to 0 for all cpus."


    overwrite: bool = False
    "Overwrite files if present."

    script: FileName | None = f"{fn}_{{species}}_{{readlength}}.sh"
    'Output Bash script. Set to "" for stdout.'


    mem: int = constants.star.indexmem
    "Memory in Gb. Set to 0 to exclude this parameter."

# star_zip_files
@app.command(group=align)
@dataclass
class star_zip_files:
help_text ="Output commands to zip star logs, outputs and fastq file lists"


    dir: FileName = "."
    "Location of input & output files"


    exec: bool = False
    "Execute command(s) instead of output."

    script: FileName | None = "unsort.sh"
    'Output commands to Bash script. Set to "" for stdout.'

# star_list_idx
@app.command(group=align)
@dataclass
class star_list_idx:
help_text ="Output list of STAR indexes found in index store(s)."

# expression counts
# featureCounts
@app.command(group=count)
@dataclass
class featureCounts:
help_text ="Run featureCounts, by default after sorting by read ID."


    bams: list[FileName], required=True
    "Input BAM file(s)."


    cpus: int =constants.cpus
    "CPUS."


    exec: bool = False
    "Execute command(s) instead of output."

    script: FileName | None = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'


    sortbyreadid: bool = True
    "Sort by read ID before running featureCounts. Use --no_sortbyreadid to skip."

    sortcpus: int int(constants.cpus / 2)
     "CPUs to use for sorting."

    sortmem: int = constants.sortmem
     "Memory to use for sorting"


    species: Literal[constants.known_species]
    "On the origin of reads."

# counts_postproc


@app.command(group=count)
@dataclass
class counts_postproc:
help_text ="Post-process featureCounts files."


    species: Literal[constants.known_species]
    "On the origin of reads."


    dir: FileName = "."
    "Location of input & output files"


    exec: bool = False
    "Execute command(s) instead of output."

    script: FileName | None = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'
# merge_counts


@app.command(group=count)
@dataclass
class merge_counts:
help_text ="Combine or extract expression values only."


    inputfiles: list[FileName], required=True
    "Input file(s)."


    dir: FileName = "."
    "Location of input & output files"


    outputdir: FileName = "."
    "Output directory."


    outputfile: FileName = constants.default_merged_counts
    "Output file."


    outputtotals: FileName = constants.default_total_counts
    "Output file for totals."


    rename_samples: FileName = constants.counts_sample_ID_file
    "Tab-delimited file with from/to pairs of sample IDs."


@app.command(group=count)
help_text ="Get sample/file IDs from featurecounts files"

@dataclass
class featurecounts_ids:


    dir: FileName | None = None
    "Dir where files of interest reside. Defaults to pwd."


    inputfiles: list[str]
    "Input file(s)."


    outputdir: str = "."
    "Output directory."


    outputfile: FileName = constants.counts_sample_ID_file
    "Output file."
# add_gene_name


@app.command(group=count)
@dataclass
class add_gene_name:
help_text ="Output expression by gene name if unique"


    species: Literal[constants.known_species]
    "On the origin of reads."

    inputfile: FileName, | None = None
    "Input file."

    outputfile: FileName, | None = constants.default_counts_by_gene_name
    "Output file."

@app.command(group=utils)
@dataclass
class rank_values:
help_text ="Add a rank column to a tab-delimited file, simply in order or based on a column."

    inputfile: FileName
    "Input file."


    outputfile: FileName | None = None
    "Output file, otherwise stdout"


    column: str | None = None
    "Column of values by which to rank. If omitted, uses the  order of appearance."


    newcolumn: str | None = None
    "Name of new column. By default, 'order' or ('rank by ' + column_name)"


    order: Literal["ascending", "descending"] = "ascending"
    "Sort order"


    sort: bool = False
    "Sort data before output"


    overwrite: bool = False
    "Overwrite files if present."


# transformcounts
@app.command(group=count)
@dataclass
class transformcounts:
help_text ="Calculate CPM or RPKM. Also output total counts."

    inputfile: FileName
    "Input file."


    method: Literal[ "CPM", "CPM-UQ", "CPM-UQ-log2", "RPKM", "RPKM-UQ", "RPKM-UQ-log2", "percentile"]
    "What to calculate."

    genenames: bool = False
    "Replace gene ID (ENGS/ENSMUS) with gene name if known and in gene metadata, otherwise use ENS ID."

    genetypes: list[str] | None = None
    'Gene type(s) to include. Set to all or "" for all', default=["protein_coding", "lncRNA"]


    log2_offset: float = 0.1
    "Offset added to values before taking log2."


    metadata: str = constants.default_gene_metadata
    "File with gene metadata. Can be original featurecounts output."


    outputfile: FileName | None = None
    "Output file. Defaults to {method}.txt"


    outputtotals: str = "total_counts"
    "Prefix for totals output files."


    overwrite: bool = False
    "Overwrite files if present."


    rescale_common_value: float = 1000
    "Target value of specified percentile after rescaling, in linear space. "


    rescale_pct: int = 75
    "Percentile of values used for rescaling. Calculated after excluding 0."
# sumcounts


@app.command(group=count)
@dataclass
class sumcounts:
help_text ="Calculate sum of counts for each sample."

    inputfile: FileName
    "Input file."


    outputfile: FileName | None = None
    "Output file. Stdout if omitted."


    overwrite: bool = False
    "Overwrite files if present."
# genebams


@app.command(group=bamfiles)
@dataclass
class genebams:
help_text ="Make BAM files with gene-specific alignments."


    bams: list[str]
    "Input BAM file(s)."

    regions: str
    "Gene regions file (bed)."


    exec: bool = False
    "Execute command(s) instead of output."


    zipbams: bool = False
    "Make a zip file of all gene bam files"

    script: str | None = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'

# rnaspades


@app.command(group=assembly)
@dataclass
class rnaspades:
help_text ="Assemble transcripts from gene-specific and unaligned reads with RNAspades."


    bams: list[str], required=True
    "Input BAM file(s)."

    cpus: int, | None = constants.cpus - 1
    "CPUs."


    exec: bool = False
    "Execute command(s) instead of output."


    genebams: list[str], | None = None,
    "Gene BAM file(s) if already created (*genes.bam)."


    mem: int = constants.rnaspades.sortmem
    "Memory in Gb."

    regions: str | None = None
    "Gene regions file (bed)."


    tempdir: str = constants.rnaspades.tempdir
    "Dir for temp files."


    zipbams: bool = False
    "Make a zip file of all gene bam files"

    script: str = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'

# abra2


@app.command(group=assembly)
@dataclass
class abra2:
help_text ="Realign indels from BAM files with abra2."

    regions: str required=True
    "Gene regions file (bed)."


    bams: list[str], | None = None,
    "BAM file(s) with all alignments (not *genes.bam)."

    cpus: int, | None = min(8, constants.cpus)
    "CPUs."


    exec: bool = False
    "Execute command(s) instead of output."


    genebams: list[str], | None = None,
    "Gene BAM file(s) if already created (*genes.bam)."
subparser.add_argument( f"--readtype: str | None = None
    "Read pairing", Literal[constants.read_types]
#Literal["paired", "single"]


    species: str | None = None, Literal[constants.known_species]
    "On the origin of reads."


    zipbams: bool = False
    "Make a zip file of all gene bam files"

    script: str | None = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'
# BAM files
# species


@app.command(group=bamfiles)
@dataclass
class species:
help_text ="Output species for every BAM file."


    bams: list[str], required=True
    "Input BAM file(s)."


    outputfile: FileName, | None = None
    "Output file. Stdout if omitted."
# readtype


@app.command(group=bamfiles)
@dataclass
class readtype:
help_text ="Output read type (paired or single) for every BAM file."


    bams: list[str], required=True
    "Input BAM file(s)."


    outputfile: FileName, | None = None
    "Output file. Stdout if omitted."
# splitbyreadgroup


@app.command(group=bamfiles)
@dataclass
class splitbyreadgroup:
help_text ="Split BAM file(s) by read group."


    bams: list[str], required=True
    "Input BAM file(s)."

    cpus: int, | None = constants.cpus
    "CPUs."

    script: str | None = f"{fn}.sh"
    'Output commands to bash script. Set to "" for stdout.'
# unsort


@app.command(group=bamfiles)
@dataclass
class unsort:
help_text ="Sort BAM by read ID."


    bams: list[str], required=True
    "Input BAM file(s)."


    exec: bool = False
    "Execute command(s) instead of output."

    script: str | None = "unsort.sh"
    'Output commands to bash script. Set to "" for stdout.'

    sortcpus: int, | None = int(constants.cpus / 2)
    "CPUs."
subparser.add_argument("--sortmem: int, | None = constants.sortmem
    "Memory.")
# junctions


@app.command(group=assembly)
@dataclass
class junctions:
help_text ="Output command to extract junctions with regtools then sum counts for identical coordinates"


    bams: list[str], required=True
    "BAM file(s)."


    mincount: int = 0
    "Minimum count"

    regions: str | None = None
    "Gene regions file (bed)."

    ext: str | None = "junctions.txt"
     "Extension for output files."

    cpus: int, | None = constants.cpus
    "CPUs."


    species: str | Literal[constants.known_species]
    "On the origin of reads."

# minimap


@app.command(group=assembly)
@dataclass
class minimap:
help_text ="Align fastas to genome"


    inputfiles: list[str]
    "Input file(s)."


    species: Literal[constants.known_species]
    "On the origin of reads."

    cpus: int, | None = constants.cpus
    "CPUs."


    exec: bool = False
    "Execute command(s) instead of output."


    format: str = "bam", Literal["sam", "bam"]
    "Output file type."


    mem: str = "3G"
    "Memory for each input batch"

    script: str | None = "minimap.sh"
    'Output commands to bash script. Set to "" for stdout.'
# other utils
# text_to_fasta


@app.command(group=utils)
@dataclass
class text_to_fasta:
help_text ="Convert tab-delimited file to fasta."

    inputfile: FileName, | None = None
    "Input file."


    minlength: int = 0
    "Minimum length"


    outputfile: FileName, | None = None
    "Output file. Stdout if omitted."

# fasta_to_text


@app.command(group=utils)
@dataclass
class fasta_to_text:
help_text ="Convert fasta file to tab-delimited text."

    inputfile: FileName, | None = None
    "Input file."


    minlength: int = 0
    "Minimum length"


    outputfile: FileName, | None = None
    "Output file. Stdout if omitted."
# join


@app.command(group=utils)
@dataclass
class join_files:
help_text ="Join two files on specified columns."


    inputfiles: list[str]
    "Input file(s)."


    columns: list[str] , | None = None, , default=["Sample_geo_accession", "experiment_alias"]
    "Column(s) to use, as one (if identical) or two comma-separated lists."


    dedup: bool = False
    "Deduplicate columns after merging."


    method: Literal["inner outer left right".split()] = "inner"
    "How to join."


    outputfile: FileName | None = None
    "Output file. Stdout if omitted."
# cat


@app.command(group=utils)
@dataclass
class concat_files:
help_text ="Concatenate files as dataframes, handle different columns"

    inputfiles: list[str]
    "Input file(s)."

    outputfile: FileName | None = None
    "Output file. Stdout if omitted."

