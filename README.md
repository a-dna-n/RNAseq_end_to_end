## What is the code in this repo?

A suite of tools that retrieves RNA-seq reads and metadata, aligns, quantifies, realigns locally etc., with just a few commands. It deals with metadata, fastq files, etc., and stores the scripts with the results.

## Motivation

1. Because alignments matter _per se_, and they're not available from Archs4, GEO recounts etc.
2. For reproducibility, transparency, consistency. > 1,200 samples  re-analyzed from scratch on modest equipment, mostly human and mouse. Manuscript in preparation.
3. But really, for the parents out there whose kids have cancer.

## Workflow

This should be a diagram, but here is a summary of the workflow:

0. Retrieve and set up reference genome & transcriptome, extract gene coordinates and other ancillary files, etc. (All versions etc. are explicit and noted.)

1. Starting from a study ID reported in an article of interest, typically a GEO series ID (GSE####), find the metadata available in GEO, SRA and ENA for samples and fastq files. Uses a combination of [pySRAdb](https://github.com/saketkc/pysradb) and direct requests.

1. Keep the original metadata files but filter, reformat and separate by species (mouse or human) and experiment type (bulk RNA-seq). Some authors have made this a royal pain.
1. Merge sample details with fastq file details; simplifies the contents, drops redundant columns, calculates the actual read length and type.
1. Choose the samples to analyze. This is one of the few manual steps, but it amounts to taking a subset of rows and columns from the last step. The result is a very simple manifest (three columns for paired reads, two for single reads). There can be multiple rows per sample, samples can be renamed,etc.
1. Make a read-length-specific index for the target species if we don't have one already.
1. Make a BED file of the target gene(s) to post-process. Technically, this is optional but it's the whole point. This command needs the species and gene name(s), but it's case-insensitive and it can optionally extend the boundaries; this distance and the species get baked into the name of the scripts that gets these boundaries.
1. Almost there: we need one more command to generate all the scripts and directories, retrieve, align etc.
    - Required:
      - the path to the manifest.
      - a directory with the STAR index to use.
    - Optionally:
      - whether to run Abra2 and/or RNAspades, which requires the file with gene boundaries and leads to small working BAM files for downstream steps0;
      - where to move the results when they're done;
      - CPUs and/or max mem for each step, where to get the fastqs (all normally based on fastq URLs, usually ENA) and how fast; everything has a default.
1. Error checking follows, based on the options (exes found, write privileges etc.). This makes temporary and destination folders, then generates separate scripts for Abra2/Spades and one wrapper Bash script that This will run everything through file cleanup, with error checking along the way, and ends with results for a group of related samples in one dir.
1. Optional: if the study authors posted any files in GEO (i.e., derived results), generate a script to retrieve them.
1. Optional: if there are NCBI counts from GEO (human only, not all studies), generate a separate script to retrieve them. Each of these scripts will download the remote files to the directory where it resides.

- Other subfunctions:


- Run Abra2 (local realignment) and RNAspades (assembly) for specific (gene) regions right after alignments.

To avoid any doubts about the provenance of the results and the methods behind them, the Bash scripts that launched most of the steps are stored with the results. They're all generated in Python.

## Status

To shorten the path to publication and make it possible for others to run this code, I'm:
- separating the functions into modules (metadata, alignment/quantification/assembly pipeline, extract specifc results from different files, making figures etc.);
- belatedly automating as many steps as possible, and
- converting all argument handling to cyclopts, like so for metadata:

<a href="docs/metadata_commands.png"><<img width=50% src="docs/metadata_commands.png"></a>

Many of these steps run in succession automatically, but heck.


## Acknowledgments

- All the diligent people who work at NCBI and ENA.
- The authors of STAR, featureCounts, Trackplot, samtools, pySRAdb, Archs4, and more.
- The study authors who posted their reads in GEO, SRA or ENA.

## Limitations

Some steps and assumptions in this pipeline are admittedly targeted adaptations, not practices I would use at work.

- For most studies of interest, read type and length were identical for all samples and fastq files were available from ENA; when there were multiple read lengths, the results were intentionally separated. For the sake of flexibility and simplicity, the same STAR index is used for all samples in a given alignment run/batch (i.e., retrieve reads, align, quantify, realign locally, assemble locally, etc.), the index is specified manually on the command line (but then recorded in a script), and there is no explicit test of the read length.
- The source of fastq files is also specified on the command line, and the URLs must be formatted according to the source, albeit simply. This is really unnecessary.
- The code that aggregates/reformats/transforms featureCounts output files involves multiple calls and doesn't get copied with the resuls.

