## What is the code in this repo?

A suite of tools that starts from an ID (GEO, SRA or ENA) and ends with alignments, expression results and assembled transcripts. Deals with metadata, fastq files, etc.

## Status

Work in progress - converting all argument handling to cyclopts.

## Motivation

Reproducibility, transparency. ~1,300 RNA-seq samples re-analyzed from scratch on modest equipment. Manuscript in preparation.

## Features (partial list)

- With two commands:
  - obtains metatadata from GEO and ENA with ([pySRAdb](https://github.com/saketkc/pysradb)) and optionally direct requests, merges and simplifies the contents (e.g. dropping uninformative or redundant columns etc.)
  - looks in GEO for supplementary files posted by study authors, and generates a Bash script to retrieve them;
  - looks in GEO for NCBI-generated counts for human experiments, and generates a Bash script to retrieve them;
  - obtains and parses the GEO sample matrix file;
  - outputs a well-formatted text file that makes it easy to choose samples to process, to know if ENA has fastqs, whether the reads are paired, and what the actual read length is.

- The user edits this file down to two or three columns: an arbitrary sample ID and the partial URL to the fastq(s).

- This simple file and one more command lead to a script that runs everything downstream (through file cleanup): retrieval, alignment, quantification, etc.
    - The key parameters for this command are the species, read length, max download speed (as a courtesy to ENA), genes of interest, and optionally network drive destination.
      - Caveat: the quantification steps involve multiple calls to the same script and are not in line with the principle that the entire procedure should be easy to follow. They should be a single Python script.

- Other subfunctions:

- Retrieve human or mouse reference genome & transcriptome, generating ancillary files e.g. gene coordinates.
- Merge featureCounts results for  groups of related samples, generate CPM and FPKM.
- Run Abra2 (local realignment) and RNAspades (assembly) for specific (gene) regions right after alignments.


## Bash from Python

To avoid any doubts about the provenance of the results and the methods behind them, the Bash scripts that launched most of the steps are stored with the results. That would be a lot of Bash scripts, but they're primarily generated in Python.

## Acknowledgments

- All the heroes at GEO, SRA and ENA.
- The authors of STAR, featureCounts, Trackplot, samtools, pySRAdb, and more.