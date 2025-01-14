## What is the code in this repo?

A suite of tools that starts from an ID (GEO, SRA or ENA) and ends with alignments, expression results and assembled transcripts. Deals with metadata, fastq files, etc.

## Motivation

Reproducibility, transparency. ~1,300 RNA-seq samples re-analyzed from scratch on modest equipment. Manuscript in preparation.

## Features (partial list)

- With two commands, obtains metatadata from both GEO and ENA with pySRAdb ([repo](https://github.com/saketkc/pysradb)), merges and simplifies them, looks in GEO for supplementary files from study authors and GEO-generated counts for human, and makes a Bash script to get them, obtains and parses the GEO sample matrix, and makes a well-formatted text file that makes it easy to choose samples to process, to know if ENA has fastqs, and what the actual read length is.
- User edits this file down to two or three columns: an arbitrary sample ID  and the partiial URL to the fastq(s).
- This simple file and one command lead to a script that runs everything downsream (through file cleanup).
- The key parameters are the species, read length, max download speed (courtesy to ENA), genes of interest, and optionally network drive destination.
- Retrieving human or mouse reference genome & transcriptome, generating ancillary files e.g. gene coordinates.
- Merging featureCounts results for  groups of related samples, generating CPM and FPKM.
- Running Abra2 and RNAspades for specific regions.

## Bash from Python

To avoid any  about the provenance of the results and the methods behind them, the Bash scripts that launched most of the steps are stored with the results. That would be a lot of Bash scripts, but they are primarily generated in Python.

## Acknowledgments

- All the heroes at GEO, SRA and ENA.
- The authors of STAR, featureCounts, Trackplot, samtools, pySRAdb, and more.