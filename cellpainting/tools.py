#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Tools
#####

*Created on Thu Jun  7 14:45 2017 by A. Pahl*

Helper Tools acting on individual data..
"""

import os
import os.path as op
from collections import Counter

import pandas as pd

from .config import ACT_PROF_PARAMETERS

ROWS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P",
        "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF"]
STRUCT = "/home/pahl/comas/share/export_data_b64.csv.gz"
KEEP = ['Compound_Id', "Batch_Id", "Producer", "Address", "Conc_uM", "Smiles", "Pure_Flag"]


def profile_sim(current, reference):
    """Calculate the similarity of two byte activity_profiles of the same length.

    Returns value between 0 .. 1"""

    ref_len = len(reference)
    assert ref_len == len(current), "Activity Profiles must have the same length to be compared."
    matching_bytes = 0
    for idx in range(ref_len):
        if current[idx] == reference[idx]:
            matching_bytes += 1
    return matching_bytes / ref_len


def format_well(well):
    """Fix well format, e.g. `A1` --> `A01`."""
    wl = len(well)
    assert wl >= 2 and wl <= 4, "well has to have 2 - 4 characters!"
    column = []
    row = []
    for pos in range(wl):
        c = well[pos]
        if c.isalpha():
            row.append(c.upper())
            continue
        row_str = "".join(row)
        assert row_str in ROWS, "row {} is not a valid row.".format(row_str)
        column.append(c)
    if len(column) < 2:
        column.insert(0, "0")  # prepend a zero
    result = row
    result.extend(column)
    return "".join(result)


def well_from_position_single(row, col):
    result = [ROWS[row - 1], "{:02d}".format(col)]
    return "".join(result)


def position_from_well_single(well):
    wl = len(well)
    column = []
    row = []
    for pos in range(wl):
        c = well[pos]
        if c.isalpha():
            row.append(c.upper())
            continue
        row_str = "".join(row)
        try:
            row_num = ROWS.index(row_str) + 1
        except ValueError:
            raise ValueError("row {} is not a valid row.".format(row_str))
        column.append(c)
    column_num = int("".join(column))
    return row_num, column_num


def find_dups(it):
    """Find duplicates in an iterable."""
    ctr = Counter(it)
    result = {}
    for c in ctr:
        if ctr[c] > 1:
            result[c] = ctr[c]
    return result


def diff(it1, it2):
    """Find the differences between two iterables"""
    s2 = set(it2)
    diff = [x for x in it1 if x not in s2]
    return diff


def print_dir(obj):
    for f in dir(obj):
        if not f.startswith("_"):
            print(f)


def create_dirs(path):
    if not op.exists(path):
        os.makedirs(path)


def check_df(df, fn):
    if df is None:
        df = pd.read_csv(fn, sep="\t")  # load default file (REFERENCES or COMAS)
    elif isinstance(df, str):
        df = pd.read_csv(df, sep="\t")
    return df


def parameters_from_act_profile_by_val(act_prof, val, parameters=ACT_PROF_PARAMETERS):
    result = []
    if not isinstance(val, str):
        val = str(val)
    for idx, act in enumerate(act_prof):
        if act == val:
            result.append(parameters[idx])
    return result
