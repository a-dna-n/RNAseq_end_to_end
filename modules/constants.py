import sys
import os
from types import SimpleNamespace
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *


class define_constants:

    #known_species = Literal["human", "mouse"] #, "monkey", "vero"]

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
        known_species = ["human", "mouse"] # , "green_monkey"] #, "monkey", "vero"]
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
        rnaspades = SimpleNamespace(
            exe = "${HOME}/bin/SPAdes-3.15.5/bin/rnaspades.py",
            tempdir = {"t3": "/media/2tb2/spades_temp"}.get(hostname, ""),
            mem = 60, #{"L490": 42, "t3": 60, "p43s": 20}.get(hostname, 0),
            sortmem = 4000 # {"L490": 4000, "t3": 4000, "p43s": 3000}.get(hostname, 0)
        )
        """
        sortcpus = 11 #{"t3": 11}.get(hostname, cpus - 1)
        bed12_columns = "chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts".split()
        gene_bams_zipfile = "gene.bams.zip"
        gene_bams_realigned_zipfile = "gene.bams.realigned.zip"
        ref = SimpleNamespace(
            dna = {
                "human": "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
                "mouse": "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
            },
            rna = {
                "human": "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz",
                "mouse": "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
            }
        )
        gtf = {
            species: os.path.join(ref_dir, os.path.basename(ref.rna[species])).replace(".gz","") for species in ref.rna.keys()
        }
        """
        minimap = SimpleNamespace(
            junctions = {
                species: os.path.basename(gtf[species]).replace("gtf", "bed") for species in gtf.keys()
            },
            genome = {
                "human": "Homo_sapiens.w5.mmi",
                "mouse": "Mus_musculus.w5.mmi",
            }
        )
        """
        species_unalias = {"Homo sapiens" :  "human", "Mus musculus": "mouse"} #, "Chlorocebus sabaeus" : "green_monkey"}

        species_tax_code = {"9606": "human", "10090" :  "mouse"} #, "60711" : "green_monkey"}
        metadata_sources = ["SRA", "GEO", "ENA"]
        file_prefix = SimpleNamespace(pysradb = "pySRAdb", geoparse = "GEOparse",  geo = "GEO",  ena = "ENA")
        self.__dict__.update(**{k: v for (k, v) in locals().items() if not(k.startswith("_")  or k.startswith("self"))})
        """
        for k, v in locals().items():
            if k.startswith("self"):
                continue
            if k.startswith("_"):
                conti   nue
            print(f"{k}: {v}")
        """

