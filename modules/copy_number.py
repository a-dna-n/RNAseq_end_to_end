#!/usr/bin/env python

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from modules.tools import *

test_libraries(["cyclopts",], exit_on_error = True)

from collections import Counter
from types import SimpleNamespace
import re
import pandas as pd
import cyclopts
from dataclasses import dataclass, KW_ONLY, field

cyc_app = cyclopts.App(help = "Functions for copy-number analyses.")
cyc_group = cyclopts.Group.create_ordered("copy_number")

nan_header = "nan"

def get_study(*, file: FileName, cna_patient_ids: list[str], unmapped: str):
    temp = pd.read_csv(file, sep="\t", header=None, dtype=object, nrows=1)
    if len(temp.columns) < 2:
        log_message("The mapping of patient to study must have two columns and be tab-delimited.", fatal=True)
    patient_to_study = pd.read_csv(file, header = None, sep="\t", usecols = temp.columns[0:2])
    (c0, c1) = patient_to_study.columns
    study_patient_ids = list(patient_to_study[c0])
    patient_to_study = patient_to_study.set_index(c0).to_dict()[c1]
    if set(cna_patient_ids) & set(study_patient_ids):
        if missing := set(cna_patient_ids) - set(study_patient_ids):
            N = len(missing)
            log_message(f"{N} patient IDs in data are not mapped to a study. They are assigned to study '{unmapped}'.", *list(missing), sep="\n\t")
            #log_message("\n\t".join([""]+list(missing)))
        if missing := set(study_patient_ids[1:]) - set(cna_patient_ids):
            N = len(missing)
            #missing = "\n\t".join(missing)
            log_message(f"{N} patients are in a study but do not have data.")
    else:
        log_message("Patient IDs in data and metadata files do not match.", fatal=True)

    new_columns = [patient_to_study.get(p, unmapped) for p in cna_patient_ids]
    #temp = "\n".join([f"\t{x}\t{y}" for (x, y)  in zip(cna_patient_ids, study)])
    temp = [f"{x}\t->\t{y}" for (x, y) in zip(cna_patient_ids, new_columns)]
    if len(temp) > 10:
        temp = temp[0:2] + ["..."] + temp[-2:]
    log_message("IDs before and after mapping:", *temp, sep ="\n\t")

    return patient_to_study

def count_cna_v1(args):
    if args.exclude:
        return count_cna_pd(args)
    verify_files_present(args.inputfile)
    if outputfile := args.outputfile:
        verify_files_absent(outputfile)
    if args.study:
        by_study = True
        verify_files_present(args.study)
    else:
        by_study = False
    cna_patient_ids = pd.read_csv(args.inputfile, header = 0, sep="\t", nrows=1).columns
    first_data_column = 1
    for c in cna_patient_ids[1:]:
        if c in args.skip:
            log_message(f"Skipping column {c}")
            first_data_column += 1
        else:
            break
    cna_patient_ids = cna_patient_ids[first_data_column:]
    log_message(f"First column of data: {cna_patient_ids[0]}")
    if by_study:
        patient_to_study = get_study(file = args.study, cna_patient_ids = cna_patient_ids, unmapped = args.unmapped)
        study = [patient_to_study.get(p, args.unmapped) for p in cna_patient_ids]
    #results = []
    unique_values_observed = set()
    bin_factor = args.bin
    with (open(outputfile, 'w') if outputfile else sys.stdout) as tempout:
        if by_study:
            print(*"gene study copy_number count", sep="\t", file=tempout)
        else:
            print(*"gene copy_number count", sep="\t", file=tempout)
            #print("gene\tcopy_number\tcount", file=tempout)
        with open(args.inputfile, encoding="utf-8") as f:
            temp = next(f).rstrip()
            for line in f:
                data_row = line.rstrip().split("\t")
                gene = data_row[0]
                data_row = data_row[first_data_column:]
                if bin_factor:
                    data_row = [x if x == "NA" else str(round(float(x) * bin_factor, 0) / bin_factor) for x in data_row]
                    #data_row = [x if x == "NA" else str(round(float(x) * bin_factor, 0) / bin_factor) for x in data_row]
                    # fix occasional -0.0
                    data_row = [x[1:] if re.match(r"^-0[\.0+]*$", x) else x for x in data_row]
                #if rounding:
                #    data_row = [x if x == "NA" else str(round(float(x), args.round)) for x in data_row]
                # r'^\-*\d+[\.\d+]*$' fails for brain_cptac_2020/data_linear_CNA.txt due to -3.81429e-05
                #data_row = [str(round(float(x) * bin_factor, 0) / bin_factor) if re.match(r'^\-*\d+[\.\d+]*$', x) else x for x in data_row]
                unique_values_observed |= set(data_row)
                if by_study:
                    counts = Counter([f"{x}\t{y}" for (x,y) in zip(study, data_row)])
                else:
                    counts = Counter(data_row)
                for k, v in counts.items():
                    print(f"{gene}\t{k}\t{v}", file=tempout)
                #results.append([gene, temp])
                #print(f"{gene}\t{repr(temp)}")
    #if by_study:
    #    if problematic := [value for value in unique_values_observed if study_separator in value]:
    #        log_message(f"values observed contain study_separator {study_separator}:" + "\n\t".join([""]+problematic), fatal=True)
    log_message("unique values observed:", *sorted(unique_values_observed), sep="\t")



@cyclopts.Parameter(name="*")
@dataclass
class _count_cna_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    study : FileName | None = None
    "File with mapping of patient to study/group."

    unmapped : str = "unk",  
    "If --study is specified, the tag assigned to any unmapped patient IDs."

    outputfile : str | None = None
    "Output file."

    bin : int = 0
    "How to bin the values, e.g. 2 for nearest 0.5 [round(x*2,0]/2"

    skip : list[str] = field(default_factory=lambda: ["Entrez_Gene_Id", "Cytoband", "Locus ID"])
    "Column(s to skip if found between gene and data. (--skip a b \"c d\""

    exclude : FileName | None = None
    "File with samples to exclude."

@cyc_app.command(group=cyc_group)
def count_cna(args: _count_cna_args) -> None:
    """Summarize quantized copy number data by gene, optionally by study/group. Reads input line by line."""
    verify_files_present(args.inputfile)
    if args.exclude:
        verify_files_present(args.exclude)
    if outputfile := args.outputfile:
        verify_files_absent(outputfile)
    if args.study:
        by_study = True
        verify_files_present(args.study)
    else:
        by_study = False
    cna_columns = pd.read_csv(args.inputfile, header = 0, sep="\t", nrows=1).columns
    if drop := set(cna_columns) & set(args.skip):
        cna_columns = [x for x in cna_columns if not x in drop]
        for c in drop:
            log_message(f"Skipping column {c}")
    if args.exclude:
        drop = set(pd.read_csv(args.exclude, header = None, sep="\t", dtype=object)[0])
        common = drop & set(cna_columns)
        log_message(f"{len(drop)} sample IDs to exclude, {len(common)} found in data")
        cna_columns = [x for x in cna_columns if not x in drop]
    data = pd.read_csv(args.inputfile, header=0, sep="\t", dtype=object, usecols=cna_columns, keep_default_na=False)
    data = data.astype(str)
    gene_col = data.columns[0]
    gene_ctr = Counter(list(data[gene_col]))
    dups = [x for x in gene_ctr.keys() if gene_ctr[x] > 1]
    if dups:
        #log_message(f"Warning: {len(dups)} genes occur more than once and be dropped:\n\t" + "\n\t".join(dups))
        log_message(f"Warning: {len(dups)} genes occur more than once and be dropped:", *dups, sep="\n\t")
        dups = set(dups)
        mask = data[gene_col].isin(dups)
        data = data[~mask].copy()
    cna_patient_ids = list(data.columns[1:])
    log_message(f"First column of data: {cna_patient_ids[0]}")
    if by_study:
        patient_to_study = get_study(file = args.study, cna_patient_ids = cna_patient_ids, unmapped = args.unmapped)
        study = [patient_to_study.get(p, args.unmapped) for p in cna_patient_ids]
        if len(set(study)) == 1:
            log_message(f"Ignore study since all samples are in {study[0]}")
            by_study = False

    unique_values_observed = set()
    bin_factor = args.bin
    with (open(outputfile, 'w') if outputfile else sys.stdout) as tempout:
        if by_study:
            print("gene\tstudy\tcopy_number\tcount", file=tempout)
        else:
            print("gene\tcopy_number\tcount", file=tempout)
        for i, data_row in data.iterrows():
            gene = data_row.values[0]
            data_row = data_row.values[1:]
            if bin_factor:
                data_row = [x if x == "NA" else str(round(float(x) * bin_factor, 0) / bin_factor) for x in data_row]
                # fix occasional -0.0
                data_row = [x[1:] if re.match(r"^-0[\.0+]*$", x) else x for x in data_row]
                # r'^\-*\d+[\.\d+]*$' fails for brain_cptac_2020/data_linear_CNA.txt due to -3.81429e-05
                #data_row = [str(round(float(x) * bin_factor, 0) / bin_factor) if re.match(r'^\-*\d+[\.\d+]*$', x) else x for x in data_row]
            #if rounding:
            #    data_row = [x if x == "NA" else str(round(float(x), args.round)) for x in data_row]
            unique_values_observed |= set(data_row)
            if by_study:
                counts = Counter([f"{x}\t{y}" for (x,y) in zip(study, data_row)])
            else:
                counts = Counter(data_row)
            for k, v in counts.items():
                print(f"{gene}\t{k}\t{v}", file=tempout)
            #results.append([gene, temp])
            #print(f"{gene}\t{repr(temp)}")
    #if by_study:
    #    if problematic := [value for value in unique_values_observed if study_separator in value]:
    #        log_message(f"values observed contain study_separator {study_separator}:" + "\n\t".join([""]+problematic), fatal=True)
    #log_message("\t".join(["unique values observed:"] + list(unique_values_observed)))
    log_message("unique values observed:", *sorted(unique_values_observed), sep="\t")


@cyclopts.Parameter(name="*")
@dataclass
class _transpose_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    outputfile : str | None = None
    "Output file."

    sort: bool = True
    "Sort output by minimum copy number column. Use --no-sort to override."


@cyc_app.command(group=cyc_group)
def transpose(args: _transpose_args):
    """Format counts with copy number values in columns."""
    verify_files_present(args.inputfile)
    if args.outputfile:
        verify_files_absent(args.outputfile)
    data = pd.read_csv(args.inputfile, header = 0, sep="\t") #, dtype={"copy_number":str, "count": int})
    data["copy_number"] = data["copy_number"].astype(str)
    unique_values_observed = list(set(data["copy_number"]))
    #if args.sort:
    #    if not args.sort in unique_values_observed:
    #        log_message(f"{args.sort} is an invalid copy number value", fatal=True)
    i = list(data.columns).index("copy_number")
    if i == 1:
        index = "gene"
    else:
        index=list(data.columns[:i])
    data = data.pivot(index=index, columns="copy_number", values="count").fillna(0).reset_index()
    new_order = list(data.columns[:i])
    numbers = [x for x in unique_values_observed if re.match(r"^-*\d+(\.\d+)*$", x)]
    temp = [[float(x), x] for x in numbers]
    temp = sorted(temp, key=lambda x: (x[0]))
    min_copy_number = temp[0][1]
    new_order += [x[1] for x in temp]
    new_order += list(set(unique_values_observed) - set(numbers))
    data = data[new_order]
    for c in unique_values_observed:
        data[c] = data[c].astype(int)
    if args.sort:
        data.sort_values(by=min_copy_number, inplace=True, ascending=False)
    data.to_csv(args.outputfile or sys.stdout, sep="\t", index=False)

@cyclopts.Parameter(name="*")
@dataclass
class _separate_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    column: str
    "Name of column to use as output file suffix"

    maxoutputfiles: int = 10
    "Maximum number of output files"

    confirm: bool = False
    "Specify this flag after checking output file names."


@cyc_app.command(group=cyc_group)
def separate(args: _separate_args):
    """Separate a file by a column, e.g. cancer type."""
    verify_files_present(args.inputfile)
    data = pd.read_csv(args.inputfile, header = 0, sep="\t", dtype=object, nrows=1)
    if args.column not in data.columns:
        log_message(f"{args.column} not found in columns of {args.inputfile}", fatal=True)
    data = pd.read_csv(args.inputfile, header = 0, sep="\t", dtype=object)
    values = set(data[args.column])
    if len(values) > args.maxoutputfiles:
        log_message(f"{len(values)} unique values in {args.column} > {args.maxoutputfiles} ", fatal=True)
    alias = {k: k.replace(" ", "_") for k in values}
    outputfile = {k: f"{args.inputfile}.{alias[k]}" for k in values}
    verify_files_absent(outputfile.values())
    for k in values:
        log_message(f"rows with {args.column} = {k} => {outputfile[k]}")
    if args.confirm:
        for i, subsetdata in data.groupby(args.column, as_index=False):
            log_message(f"Writing to {outputfile[i]}")
            subsetdata.to_csv(outputfile[i], sep="\t", index=False)
    else:
        log_message(f"Repeat with --confirm if file names are OK.")


def _exclude_nan_column(*, data, keep_total: bool=False) -> pd.DataFrame:
    if not nan_header in data.columns:
        log_message(f"{nan_header} column not found")
        return data
    data.reset_index(inplace=True)
    if not "gene" in data.columns:
        log_message(f"Unknown data", fatal=True)
    if data.columns[1] == "study":
        data.set_index(["gene", "study"], inplace=True)
    else:
        data.set_index("gene", inplace=True)
    data["total"] = data.sum(axis=1)
    mask = (data["total"] ==  data[nan_header]) & (data[nan_header] > 0)
    all_nans = list(data[mask].index) #["gene"])
    data = data[~mask].copy()
    if all_nans:
        all_nans = sorted(all_nans)
        log_message(f"Dropping {len(all_nans)} genes with 100% nan:", *all_nans, sep = "\n")
    else:
        log_message("No genes with 100% nans")
    if keep_total:
        data["total"] = data["total"] - data[nan_header]
    data.drop(columns = nan_header, inplace=True)
    return data #.reset_index()

@cyc_app.command(group=cyc_group)
def drop_nans(*, inputfile: FileName, outputfile: FileName | None = None):
    verify_files_present(args.inputfile)
    if args.outputfile:
        verify_files_absent(args.outputfile)
    data = pd.read_csv(args.inputfile, header = 0, sep="\t")
    data = _exclude_nan_column(data=data, keep_total=False)
    data.to_csv(args.outputfile if args.outputfile else sys.stdout, sep="\t")


@cyclopts.Parameter(name="*")
@dataclass
class _fraction_args:

    _: KW_ONLY

    inputfile: FileName
    "Input file."

    outputfile : str | None = None
    "Output file."

    round: int = 2
    "Round to this many decimal places."

    exclude_nan: bool = True
    "Exclude nan counts, drop but report any genes with 0 counts remaining."

@cyc_app.command(group=cyc_group)
def fractions(args: _fraction_args):
    verify_files_present(args.inputfile)
    if outputfile := args.outputfile:
        verify_files_absent(outputfile)

    data = pd.read_csv(args.inputfile, header = 0, sep="\t")
    if exclude_nan:
        data = _exclude_nan_column(data=data, keep_total=True)
    data_columns = list(data.columns)[:-1]
    #data["total"] = data.sum(axis=1) #.reset_index()
    for c in data_columns:
        data[c] /= data["total"]
        data[c] *= 100
        data[c] = data[c].round( args.round)
    data.sort_values(by=data_columns[0], ascending=False).to_csv(outputfile if outputfile else sys.stdout, sep="\t")

"""
def counts_to_fractions(args):
    verify_files_present(args.inputfile)
    if outputfile := args.outputfile:
        verify_files_absent(outputfile)
    data = pd.read_csv(args.inputfile, sep="\t", header=0, index_col=0)
    columns = data.columns
    if len(set(data.dtypes)) == 1:
        dtype_for_total = data.dtypes[0]
    else:
        dtype_for_total = None
    data["total"] = data.apply(lambda x: sum(x), axis=1).astype(int)
    for x in data.columns[:-1]:
        data[x] /= data["total"]
        if args.pct:
            data[x] *= 100
    data = data.round(args.round)
    if dtype_for_total:
        data["total"] = data["total"].astype(dtype_for_total)
    data.to_csv(outputfile or sys.stdout, sep="\t") #, index=False)

def add_fractions(args):
    verify_files_present(args.inputfile)
    if outputfile := args.outputfile:
        verify_files_absent(outputfile)
    data = pd.read_csv(args.inputfile, sep="\t", header=0, index_col=0)
    number_columns = data.columns
    if len(set(data.dtypes)) == 1:
        dtype_for_total = data.dtypes[0]
    else:
        dtype_for_total = None
    data["total"] = data.apply(lambda x: sum(x), axis=1).astype(int)
    for col in number_columns:
        pct_col = f"%({col})"
        data[pct_col] = round(100 * data[col] / data["total"], args.round)
    if  args.drop_COUNTS:
        data.drop(columns = number_columns, inplace=True)
    if args.sortby:
        ok = False
        col =  args.sortby.replace('"',"")
        if col.endswith("-"):
            col = col[:-1]
            _ascending = False
        else:
            _ascending = True
        if col in data.columns:
            pass
            ok = True
        elif col.isdigit() and int(col) >= 0 and int(col) < len(data.columns):
            col = data.columns[int(col)]
        else:
            log_message(f"{col} not found in data", fatal=True)
        log_message(f"sorting by {col}")
        data = data.sort_values(by=col, ascending=_ascending)
    data.to_csv(outputfile or sys.stdout, sep="\t") #, index=False)

def subset_by_row_ids(args):
    verify_files_present([args.inputfile, args.ids])
    if args.outputfile:
        verify_files_absent(args.outputfile)

    #data = pd.read_csv(args.inputfile, sep="\t", header=0, dtype=object)
    id_list = set(pd.read_csv(args.ids, sep="\t", header=None, dtype=object)[0])
    data = pd.read_csv(args.inputfile, sep="\t", header=0, dtype=object)
    mask = data[data.columns[0]].isin(id_list)
    data[mask].to_csv(args.outputfile or sys.stdout, sep="\t", index=False)

brain_cptac_2020/data_linear_CNA.txt.counts.2
-5.5	5.0	-15.0	-0.0	-4.5	-13.5	-3.5	-17.5	3.5	1.5	NA	0.0	-1.5	-6.5	4.0	6.0	-2.0	5.5	-14.0	0.5	-1.0	-18.0	-3.0	-17.0	-12.5	2.5	-10.5	-4.0	-16.5	-0.5	-5.0	1.0	2.0	4.5	3.0	-2.5
openpedcan_v15/data_CNA.txt.counts
0	-2	2	-1	1
openpedcan_v15/data_linear_CNA.txt.counts
50	4	129	10	93	72	80	84	89	105	171	142	18	121	622	125	14	138	68	127	248	144	145	21	33	32	106	59	149	246	128	4.5	95	19	54	122	140	0	35	5	28	51	30	57	63	118	53	188	96	58	152	26	112	17	76	37	2	48	16	39	119	364	85	213	20	230	117	116	38	114	238	67	203	34	43	60	240	522	82	132	599	131	147	29	75	45	103	81	329	86	222	197	3	74	56	109	55	71	42	65	195	225	61	44	130	5.5	49	24	66	77	647	25	200	216	62	36	27	181	97	13	22	40	92	113	23	8	31	52	47	64	164	69	83	102	1	70	12	46	137	3.5	6	15	170	115	284	101	244	99	260	11	100	41	73	7	9
pptc/data_CNA.txt.counts
0	-2	2	-1
x01_fy16_nbl_maris/data_CNA.txt.counts
0	1	-2	-1	2
x01_fy16_nbl_maris/data_linear_CNA.txt.counts
18	659	75	29	70	6	104	28	115	131	147	622	67	101	23	14	49	145	10	164	16	122	188	3	68	36	232	11	2	12	4	42	246	37	24	76	26	66	50	32	31	30	99	25	43	13	203	65	462	165	220	225	89	20	222	69	137	71	196	148	39	5	19	7	1	248	38	106	15	126	59	130	56	57	51	112	53	1008	8	55	45	33	284	41	774	34	60	256	230	504	58	72	132	110	81	46	21	17	0	9	22	54	27	339	329	647	35

values observed:
-2
0
1
-1
2

"""
if __name__ == "__main__":
    cyc_app()


