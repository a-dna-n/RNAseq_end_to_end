import sys
import os
import argparse
from types import SimpleNamespace
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.constants import *

# constants = define_constants()

class ArgDef:

    def __init__(
        self, *, desc: str, subfunctions: bool = True, allow_unknown_args: bool = False
    ):
        self.unknownargs = None
        self.args = None
        self.argparser = argparse.ArgumentParser(
            description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subfunctions = subfunctions
        if subfunctions:
            self.subparsers = self.argparser.add_subparsers(
                description=" ", metavar=""
            )
        else:
            self.subparsers = None
        self.subparser = {}
        self.allow_unknown_args = allow_unknown_args
        

    def add_subfunction(self, *, newfn: str, help_text: str, args: list|None = None):
        if self.subparsers is None:
            sys.exit("self.subparsers not initialized")
        subparser = self.subparsers.add_parser(
            newfn, help=help_text, formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        subparser.set_defaults(func=newfn)
        self.subparser[newfn] = subparser
        return subparser

    @staticmethod
    def convert_args_to_upper_case(data: dict | argparse.Namespace):
        """
        Converts argparse parameter names to uppercase, e.g. the value of --inputfile is accessed as args.INPUTFILE.
        Nested dictionaries not tested.

        Args:
            data: a dictionary or namespaces (not nested)

        Returns:
            A simplenamespace with dot notation.
        """
        if isinstance(data, argparse.Namespace):
            return SimpleNamespace(**{k.upper(): v for k, v in data.__dict__.items()})
        else:
            return SimpleNamespace(**{k.upper(): v for k, v in data.items()})

    def handle_args(self, assume_help: bool = True):
        if len(sys.argv) == 1 and assume_help:
            self.argparser.print_help()
            sys.exit()
        if self.allow_unknown_args:
            # knownargs, unknownargs = self.argparser.parse_known_args()
            args, self.unknownargs = self.argparser.parse_known_args()
        else:
            args = self.argparser.parse_args()
        if self.subfunctions:
            self.func = args.func
        else:
            self.func = None
        self.args = self.convert_args_to_upper_case(args)


class ArgDef:

    def __init__(
        self, *, desc: str, subfunctions: bool = True, allow_unknown_args: bool = False
    ):
        self.unknownargs = None
        self.args = None
        self.argparser = argparse.ArgumentParser(
            description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subfunctions = subfunctions
        if subfunctions:
            self.subparsers = self.argparser.add_subparsers(
                description=" ", metavar=""
            )
        else:
            self.subparsers = None
        self.subparser = {}
        self.allow_unknown_args = allow_unknown_args
        

    def add_subfunction(self, *, newfn: str, help_text: str, args: list|None = None):
        if self.subparsers is None:
            sys.exit("self.subparsers not initialized")
        subparser = self.subparsers.add_parser(
            newfn, help=help_text, formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        subparser.set_defaults(func=newfn)
        self.subparser[newfn] = subparser
        return subparser

    @staticmethod
    def convert_args_to_upper_case(data: dict | argparse.Namespace):
        """
        Converts argparse parameter names to uppercase, e.g. the value of --inputfile is accessed as args.INPUTFILE.
        Nested dictionaries not tested.

        Args:
            data: a dictionary or namespaces (not nested)

        Returns:
            A simplenamespace with dot notation.
        """
        if isinstance(data, argparse.Namespace):
            return SimpleNamespace(**{k.upper(): v for k, v in data.__dict__.items()})
        else:
            return SimpleNamespace(**{k.upper(): v for k, v in data.items()})

    def handle_args(self, assume_help: bool = True):
        if len(sys.argv) == 1 and assume_help:
            self.argparser.print_help()
            sys.exit()
        if self.allow_unknown_args:
            # knownargs, unknownargs = self.argparser.parse_known_args()
            args, self.unknownargs = self.argparser.parse_known_args()
        else:
            args = self.argparser.parse_args()
        if self.subfunctions:
            self.func = args.func
        else:
            self.func = None
        self.args = self.convert_args_to_upper_case(args)

def define_args(constants):

    args = ArgDef(
        desc="Functions for RNA-seq data.", subfunctions=True, allow_unknown_args=False
    )

    fn = "metadata"
    help_text ="Find metadata and files available from GEO, SRA and ENA."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--studydir",
        "-d",
        required=True,
        help="Folder where metadata files will be stored, with each source in its own subfolder.",
    )
    subparser.add_argument(
        "--id", "-i", type=str, required=True, help="GEO, SRA or ENA identifier. Ex. PRJNA743892 GSE162198 SRP294329",
    )
    """
    subparser.add_argument(
        "--sources",
        "-s",
        type=str,
        required=False,
        help="Source(s) of metadata to query",
        choices = constants.metadata_sources, #= ["SRA", "GEO", "ENA"],
        default = constants.metadata_sources, #["ENA", "SRA", "GEO"],
        nargs = "+"
    )
    """
    subparser.add_argument(
        "--assay",
        "-a",
        type=str,
        required=False,
        help="Data type of interest.",
        default = ["RNA-Seq"],
        nargs = "+"
    )
    subparser.add_argument(
        "--species",
        "-S",
        type=str,
        required=False,
        help="Species of interest.",
        default = constants.known_species,
        nargs = "+"
    )
    
    # metadata
    fn = "pySRAdb_get_metadata"
    help_text ="Get metadata with pysradb"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--ids", "-i", type=str, required=True, help="Identifier(s).", nargs="+"
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file for single ID. Default is {ID}.{ext}",
    )
    subparser.add_argument(
        "--ext",
        "-e",
        type=str,
        required=False,
        help="Output file suffix",
        default="_pysradb.txt",
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    # subparser.add_argument("--source", "-m", type=str, required=False, help="Source of metadata.", default="SRA", choices=["ENA", "SRA", "GEO"])

    # geo
    fn = "geo"
    help_text ="Get metadata from GEO"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--ids", "-i", type=str, required=True, help="Identifier(s).", nargs="+"
    )
    subparser.add_argument(
        "--ext",
        "-O",
        type=str,
        required=False,
        help="Output file suffix",
        default="_geo.txt",
    )

    # dedup_cols
    fn = "dedup_cols"
    help_text ="Remove duplicate columns"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    # enafastqs
    fn = "enafastqs"
    help_text ="Get a ist of ENA fastqs for ID(s) like PRJNA627881"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--ids", "-i", type=str, required=True, help="Identifier(s).", nargs="+"
    )
    subparser.add_argument(
        "--keeptempfiles",
        "-k",
        required=False,
        help="Keep temp output files.",
        default=False,
        action="store_true",
    )
    subparser.add_argument(
        "--ext",
        "-O",
        type=str,
        required=False,
        help="Output file extension",
        default="ena.txt",
    )

    # genecoords
    fn = "get_genecoords"
    help_text ="Get coordinates for a list of genes, in BED format, or output the script to do this."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--genes",
        "-g",
        type=str,
        required=True,
        help="Gene(s) in a file or as a space-separated list",
        nargs="+",
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--extend",
        "-x",
        type=int,
        required=False,
        help="Extension of gene boundaries, in nucleotides. Maximum not checked against chromosome length.",
        default=0,
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file for gene coordinates. Stdout if omitted.",
        default = "gene_coords_{species}_ext_{extend}.bed"
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        "--script", "-S", type=str, required=False, #default=f"{fn}.sh",
        help='Output command to this Bash script instead of executing.'
    )

    # gene_metadata
    fn = "gene_metadata"
    help_text ="Extract metadata columns from featureCounts output files."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=True, help="Input file."
    )
    subparser.add_argument(
        "--outputdir",
        "-O",
        type=str,
        required=False,
        help="Output directory.",
        default=".",
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file for gene metadata.",
        default=constants.default_gene_metadata,
    )

    # pre-STAR

    fn = "genomeref"
    help_text ="Generate script to fetch genome (fasta) and transcriptome (gtf) files."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--outputdir",
        "-O",
        type=str,
        required=False,
        help="Output directory.",
        default=constants.ref_dir
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}_{{species}}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )

    fn = "gtf_to_coords"
    help_text ="Extract gene coordinates from GTF."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file. Required unless --species is specified."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file.  Required unless --species is specified.",
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=False,
        help="On the origin of reads. Assumes default input & output files.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        "--sort",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Sort by chr, start & stop. Use --no_sort to skip.",
    )
    subparser.add_argument(
        "--genes", "-g", type=str, required=False, help="Gene type(s) to keep, space-separated. Use \"\" for all.", default=["protein_coding"], nargs="+"
    )

    fn = "ens_to_gene"
    help_text ="Extract unique gene-name/gene-ID pairs from GTF."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file. Required unless --species is specified."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file.  Required unless --species is specified.",
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=False,
        help="On the origin of reads. Assumes default input & output files.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )

    # STAR alignments etc.

    fn = "star"
    help_text ="Output commands to fetch fastqs then run STAR & featureCounts. Input=label tab fastq1 (tab fastq2)"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--index",
        "-x",
        type=str,
        required=True,
        help="Genome index. Path is ignored.",
        choices=list(constants.star_indexes),
    )
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=True, help="Input file."
    )
    subparser.add_argument(
        "--abra2",
        type=str,
        required=False,
        help="Coords file to use for gene bams then abra2.",
    )
    subparser.add_argument(
        "--additional",
        "-a",
        type=str,
        required=False,
        help="Additional options, in double quotes, added verbatim.",
        default="",
    )
    subparser.add_argument(
        "--addreadgroup",
        required=False,
        help="Use sample ID as read group.",
        action="store_true",
    )
    subparser.add_argument(
        "--bamdest",
        "-D",
        type=str,
        required=False,
        help="Destination for BAM files after completion, e.g. on mounted drive",
    )
    subparser.add_argument(
        "--counts",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run featurecounts after alignments. Use --no_counts to skip.",
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUs."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--fastqsource",
        "-f",
        type=str,
        required=False,
        help="Source of fastq files.",
        default="ENA",
        choices=["ENA", "SRA"],
    )
    subparser.add_argument(
        "--outputdir",
        "-O",
        type=str,
        required=False,
        help="Output directory.",
        default=".",
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        f"--readtype",
        "-r",
        type=str,
        required=False,
        help="Read pairing",
        choices=constants.read_types,
    )
    """
    subparser.add_argument(
        "--rnaspades",
        type=str,
        required=False,
        help="Coords file to use for gene bams then rnaspades.",
    )
    """
    subparser.add_argument(
        "--sortcpus", "-C", type=int, required=False, default=int(constants.cpus / 2),
        help="CPUs to use for sorting"
    )
    subparser.add_argument("--sortmem", "-M", type=int, required=False, default=2000, help="Memory to use for sorting")
    subparser.add_argument(
        "--transferspeed",
        "-t",
        type=int,
        required=False,
        help="Aspera transfer speed.",
        default=500,
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )

    fn = "star_clear_mem"
    help_text ="Output command to remove star genome from mem"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--index", "-x", type=str, required=True, help="Genome index. Path is ignored."
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )

    fn = "star_make_idx"
    help_text ="Output command to make an index for a given species and read length."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--readlength",
        "-r",
        type=int,
        required=True,
        help="Read length minus 1, probably 49, 75, 100 or  150. See Star Manual"
    )
    subparser.add_argument(
        "--index",
        "-i",
        type=str,
        required=False,
        help="Name of index, i.e., base output directory.",
        default=os.path.join(constants.star.base_dir, "{species}.{dna_build}_{rna_build}.{readlength}")
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUs. Set to 0 for all cpus."
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}_{{species}}_{{readlength}}.sh",
        help='Output Bash script. Set to "" for stdout.',
    )
    subparser.add_argument(
        "--mem",
        "-m",
        type=int,
        required=False,
        help="Memory in Gb. Set to 0 to exclude this parameter.",
        default=constants.star.indexmem,
    )

    # star_zip_files
    fn = "star_zip_files"
    help_text ="Output commands to zip star logs, outputs and fastq file lists"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--dir",
        "-d",
        type=str,
        required=False,
        help="Location of input & output files",
        default=".",
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default="unsort.sh",
        help='Output commands to Bash script. Set to "" for stdout.',
    )

    # star_list_idx
    fn = "star_list_idx"
    help_text ="Output list of STAR indexes found in index store(s)."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)

    # expression counts

    # featureCounts
    fn = "featureCounts"
    help_text ="Run featureCounts, by default after sorting by read ID."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUS."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )
    subparser.add_argument(
        "--sortbyreadid",
        "-r",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Sort by read ID before running featureCounts. Use --no_sortbyreadid to skip.",
    )
    subparser.add_argument(
        "--sortcpus", "-C", type=int, required=False, default=int(constants.cpus / 2),
        help= "CPUs to use for sorting."
    )
    subparser.add_argument(
        "--sortmem", "-M", type=int, required=False, default=constants.sortmem,
        help = "Memory to use for sorting"
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=False,
        help="On the origin of reads.",
        choices=constants.known_species,
    )

    # counts_postproc
    fn = "counts_postproc"
    help_text ="Post-process featureCounts files."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species,
    )
    subparser.add_argument(
        "--dir",
        "-d",
        type=str,
        required=False,
        help="Location of input & output files",
        default=".",
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )

    # merge_counts
    fn = "merge_counts"
    help_text ="Combine or extract expression values only."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfiles", "-i", type=str, required=True, help="Input file(s).", nargs="+"
    )
    subparser.add_argument(
        "--dir",
        "-d",
        type=str,
        required=False,
        help="Location of input & output files",
        default=".",
    )
    subparser.add_argument(
        "--outputdir",
        "-O",
        type=str,
        required=False,
        help="Output directory.",
        default=".",
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file.",
        default=constants.default_merged_counts,
    )
    subparser.add_argument(
        "--outputtotals",
        "-t",
        type=str,
        required=False,
        help="Output file for totals.",
        default=constants.default_total_counts,
    )
    subparser.add_argument(
        "--rename_samples",
        "-r",
        type=str,
        required=False,
        help="Tab-delimited file with from/to pairs of sample IDs.",
        default=constants.counts_sample_ID_file,
    )
    # subparser.add_argument("--ext", "-O", type=str, required=False, default= constants.default_merged_counts)

    # featurecounts_ids
    fn = "featurecounts_ids"
    help_text ="Get sample/file IDs from featurecounts files"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    # subparser.add_argument("--dir", "-d", type=str, required=False, help="Dir where files of interest", default=".")
    subparser.add_argument(
        "--inputfiles", "-i", type=str, required=True, help="Input file(s).", nargs="+"
    )
    subparser.add_argument(
        "--outputdir",
        "-O",
        type=str,
        required=False,
        help="Output directory.",
        default=".",
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file.",
        default=constants.counts_sample_ID_file,
    )

    # add_gene_name
    fn = "add_gene_name"
    help_text ="Output expression by gene name if unique"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species
    )
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        default=constants.default_counts_by_gene_name,
        help="Output file."
    )

    fn = "rank_values"
    help_text ="Add a rank column to a tab-delimited file, simply in order or based on a column."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=True, help="Input file."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file, otherwise stdout",
    )
    subparser.add_argument(
        "--column",
        "-c",
        type=str,
        required=False,
        help="Column of values by which to rank. If omitted, uses the  order of appearance.",
    )
    subparser.add_argument(
        "--newcolumn",
        "-n",
        type=str,
        required=False,
        help="Name of new column. By default, 'order' or 'rank by ' + column_name",
    )
    subparser.add_argument(
        "--order",
        "-O",
        help="Sort order",
        type=str,
        required=False,
        default = "ascending",
        choices = ["ascending", "descending"],
    )
    subparser.add_argument(
        "--sort",
        "-s",
        required=False,
        help="Sort data before output",
        action="store_true",
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )

    # transformcounts
    fn = "transformcounts"
    help_text ="Calculate CPM or RPKM. Also output total counts."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=True, help="Input file."
    )
    subparser.add_argument(
        "--method",
        "-M",
        type=str,
        required=True,
        help="What to calculate.",
        choices=[
            "CPM",
            "CPM-UQ",
            "CPM-UQ-log2",
            "RPKM",
            "RPKM-UQ",
            "RPKM-UQ-log2",
            "percentile",
        ],
    )
    subparser.add_argument(
        "--genenames",
        "-n",
        required=False,
        help="Replace gene ID (ENGS/ENSMUS) with gene name if known and in gene metadata, else use ENS ID.",
        action="store_true",
    )
    subparser.add_argument(
        "--genetypes",
        "-g",
        type=str,
        required=False,
        help='Gene type(s) to include. Set to all or "" for all',
        nargs="+",
        default=["protein_coding", "lncRNA"],
    )
    subparser.add_argument(
        "--log2_offset",
        "-l",
        type=float,
        required=False,
        help="Offset added to values before taking log2.",
        default=0.1,
    )
    subparser.add_argument(
        "--metadata",
        "-m",
        type=str,
        required=False,
        help="File with gene metadata. Can be original featurecounts output.",
        default=constants.default_gene_metadata,
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Defaults to {method}.txt",
    )
    subparser.add_argument(
        "--outputtotals",
        "-t",
        type=str,
        required=False,
        help="Prefix for totals output files.",
        default="total_counts",
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )
    subparser.add_argument(
        "--rescale_common_value",
        "-r",
        type=float,
        required=False,
        help="Target value of specified percentile after rescaling, in linear space.",
        default=1000,
    )
    subparser.add_argument(
        "--rescale_pct",
        "-p",
        type=int,
        required=False,
        help="Percentile to use for rescaling. Calculated after excluding 0.",
        default=75,
    )

    # sumcounts
    fn = "sumcounts"
    help_text ="Calculate sum of counts for each sample."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=True, help="Input file."
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )
    subparser.add_argument(
        "--overwrite",
        "-w",
        required=False,
        help="Overwrite files if present.",
        action="store_true",
    )

    # genebams
    fn = "genebams"
    help_text ="Make BAM files with gene-specific alignments."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--regions", "-r", type=str, required=True, help="Gene regions file (bed)."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--zipbams",
        "-z",
        required=False,
        help="Make a zip file of all gene bam files",
        action="store_true",
    )
    subparser.add_argument(
        "--script", "-S", type=str, required=False, default=f"{fn}.sh", 
        help='Output commands to bash script. Set to "" for stdout.',
    )
    """
    # rnaspades
    fn = "rnaspades"
    help_text ="Assemble transcripts from gene-specific and unaligned reads with RNAspades."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus - 1, help="CPUs."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--genebams",
        "-g",
        type=str,
        required=False,
        help="Gene BAM file(s) if already created (*genes.bam).",
        nargs="+",
    )
    subparser.add_argument(
        "--mem",
        "-m",
        type=int,
        required=False,
        help="Memory in Gb.",
        default=constants.rnaspades.sortmem,
    )
    subparser.add_argument(
        "--regions", "-r", type=str, required=False, help="Gene regions file (bed)."
    )
    subparser.add_argument(
        "--tempdir",
        "-t",
        type=str,
        required=False,
        help="Dir for temp files.",
        default=constants.rnaspades.tempdir,
    )
    subparser.add_argument(
        "--zipbams",
        "-z",
        required=False,
        help="Make a zip file of all gene bam files",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )
    """
    # abra2
    fn = "abra2"
    help_text ="Realign indels from BAM files with abra2."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--regions", "-r", type=str, required=True, help="Gene regions file (bed)."
    )
    subparser.add_argument(
        "--bams",
        "-b",
        type=str,
        required=False,
        help="BAM file(s) with all alignments (not *genes.bam).",
        nargs="+",
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=min(8, constants.cpus), help="CPUs."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--genebams",
        "-g",
        type=str,
        required=False,
        help="Gene BAM file(s) if already created (*genes.bam).",
        nargs="+",
    )
    subparser.add_argument(
        f"--readtype",
        "-t",
        type=str,
        required=False,
        help="Read pairing",
        choices=constants.read_types, #Literal["paired", "single"],
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=False,
        help="On the origin of reads.",
        choices=constants.known_species
    )
    subparser.add_argument(
        "--zipbams",
        "-z",
        required=False,
        help="Make a zip file of all gene bam files",
        action="store_true",
    )
    subparser.add_argument(
        "--script", "-S", type=str, required=False, default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.'
    )

    # BAM files

    # species
    fn = "bam_species"
    help_text ="Output species for every BAM file."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    # readtype
    fn = "readtype"
    help_text ="Output read type (paired or single) for every BAM file."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    # splitbyreadgroup
    fn = "splitbyreadgroup"
    help_text ="Split BAM file(s) by read group."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUs."
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )

    # unsort
    fn = "unsort"
    help_text ="Sort BAM by read ID."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="Input BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--script",
        "-S",
        type=str,
        required=False,
        default="unsort.sh",
        help='Output commands to bash script. Set to "" for stdout.',
    )
    subparser.add_argument(
        "--sortcpus", "-C", type=int, required=False, default=int(constants.cpus / 2), help="CPUs."
    )
    subparser.add_argument("--sortmem", "-M", type=int, required=False, default=constants.sortmem, help="Memory.")

    # junctions
    fn = "junctions"
    help_text ="Output command to extract junctions with regtools then sum counts for identical coordinates"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--bams", "-b", type=str, required=True, help="BAM file(s).", nargs="+"
    )
    subparser.add_argument(
        "--mincount", "-m", type=int, required=False, help="Minimum count", default=0
    )
    subparser.add_argument(
        "--regions", "-r", type=str, required=False, help="Gene regions file (bed)."
    )
    subparser.add_argument(
        "--ext", "-O", type=str, required=False, default="junctions.txt", help = "Extension for output files."
    )
    """
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUs."
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=False,
        help="On the origin of reads.",
        choices=constants.known_species
    )
    """
    """
    # minimap
    fn = "minimap"
    help_text ="Align fastas to genome"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfiles", "-i", type=str, required=True, help="Input file(s).", nargs="+"
    )
    subparser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="On the origin of reads.",
        choices=constants.known_species
    )
    subparser.add_argument(
        "--cpus", "-c", type=int, required=False, default=constants.cpus, help="CPUs."
    )
    subparser.add_argument(
        "--exec",
        "-X",
        required=False,
        help="Execute command(s) instead of output.",
        action="store_true",
    )
    subparser.add_argument(
        "--format",
        "-f",
        type=str,
        required=False,
        help="Output file type.",
        default="bam",
        choices=["sam", "bam"],
    )
    subparser.add_argument(
        "--mem",
        "-m",
        type=str,
        required=False,
        help="Memory for each input batch",
        default="3G",
    )
    subparser.add_argument(
        "--script", "-S", type=str, required=False, default=f"{fn}.sh",
        help='Output commands to bash script. Set to "" for stdout.'
    )

    # other utils

    # text_to_fasta
    fn = "text_to_fasta"
    help_text ="Convert tab-delimited file to fasta."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file."
    )
    subparser.add_argument(
        "--minlength", "-m", type=int, required=False, help="Minimum length", default=0
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )
    """

    # fasta_to_text
    fn = "fasta_to_text"
    help_text ="Convert fasta file to tab-delimited text."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfile", "-i", type=str, required=False, help="Input file."
    )
    subparser.add_argument(
        "--minlength", "-m", type=int, required=False, help="Minimum length", default=0
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    # join
    fn = "join_files"
    help_text ="Join two files on specified columns."
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfiles", "-i", type=str, required=True, help="Input file(s).", nargs="+"
    )
    subparser.add_argument(
        "--columns",
        "-c",
        type=str,
        required=False,
        help="Column(s) to use, as one (if identical) or two comma-separated lists.",
        nargs="+",
        default=["Sample_geo_accession", "experiment_alias"],
    )
    subparser.add_argument(
        "--dedup",
        "-d",
        required=False,
        help="Deduplicate columns after merging.",
        action="store_true",
    )
    subparser.add_argument(
        "--method",
        "-M",
        type=str,
        required=False,
        help="How to join.",
        choices="inner outer left right".split(),
        default="inner",
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    # cat
    fn = "concat_files"
    help_text ="Concatenate files as dataframes, handle different columns"
    subparser = args.add_subfunction(newfn=fn, help_text=help_text)
    subparser.add_argument(
        "--inputfiles", "-i", type=str, required=True, help="Input file(s).", nargs="+"
    )
    subparser.add_argument(
        "--outputfile",
        "-o",
        type=str,
        required=False,
        help="Output file. Stdout if omitted.",
    )

    return args
