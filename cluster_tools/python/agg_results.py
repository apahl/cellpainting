#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############################
Aggregate CellProfiler Results
##############################

*Created on Thu Aug 22 12:00 2017 by A. Pahl*"""


import argparse
import sys
import os.path as op
# from collections import Counter

import pandas as pd


def usage():
    print("""Aggregate CellProfiler Results by Mean or Median.
    Writes out a Results.tsv file.
    agg_results.py -h for more help.""")
    sys.exit(1)


def flush_print(txt):
    txt = txt + "\r"
    print(txt)
    sys.stdout.flush()


def aggregate(input_dir, agg_type="median", sep="\t"):
    df_list = []
    keep_image = ["ImageNumber", "Metadata_Well", "Metadata_Plate", "Metadata_Site", "Count_Cells"]
    if sep == ",":
        f_ext = "csv"
    else:
        f_ext = "txt"
    for idx in range(96):
        im_start = idx * 36 + 1
        flush_print("* Slice {:2d}: {} - {}...".format(idx + 1, im_start, im_start + 35))
        df_slice = pd.read_csv("{}/{}/Image.{}".format(input_dir, im_start, f_ext), sep=sep)
        df_slice = df_slice[keep_image]
        for ch in ["Cells", "Nuclei", "Cytoplasm"]:
            df = pd.read_csv("{}/{}/{}.{}".format(input_dir, im_start, ch, f_ext), sep=sep)
            keys = list(df.keys())
            keys.pop(keys.index("ImageNumber"))
            keys.pop(keys.index("ObjectNumber"))
            if agg_type == "median":
                df = df.groupby("ImageNumber")[keys].median()
                df = df.add_prefix("Median_{}_".format(ch)).reset_index()
            else:
                df = df.groupby("ImageNumber")[keys].mean()
                df = df.add_prefix("Mean_{}_".format(ch)).reset_index()
            df_slice = pd.merge(df_slice, df, on="ImageNumber", how="left")
        df_list.append(df_slice)
    df_plate = pd.concat(df_list)
    nrows = df_plate.shape[0]
    if nrows != 3456:
        print("# unexpected number of rows:", nrows)
    df_plate.to_csv("{}/Results.tsv".format(input_dir), sep="\t")


if __name__ == "__main__":
    # file_to_search file_w_smiles output_dir job_index
    parser = argparse.ArgumentParser(
        descripiton="Aggregate CellProfiler Results by Mean or Median.\nWrites out a Results.tsv file.")
    parser.add_argument("input_dir", help="The directory with the CellProfiler results.")
    parser.add_argument("-s", "--sep", help="column separator (default is tab).")
    parser.add_argument("-t", "--type", help="aggregation type, mean or median.",
                        choices=["mean", "median"])
    # parser.add_argument("", help="")

    args = parser.parse_args()
    sep = "\t"
    if args.sep is not None:
        sep = args.sep
    if args.type is None:
        args.type = "median"
    print("* aggregation type set to", args.type)
    err = False
    reason = ""
    if not op.isdir(args.input_dir):
        args.input_dir = args.input_dir + "_output"
        if not op.isdir(args.input_dir):
            err = True
            reason = "Input dir {} does not exist.".format(args.input_dir)

    if err:
        print(reason)
        usage()

    aggregate(args.input_dir, agg_type=args.type, sep=sep)
