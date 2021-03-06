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


try:
    from .nim_ext import profile_sim
    print("- loaded Nim extension.")
except ImportError:
    def profile_sim(current, reference, ignore_direction=False, mask=""):
        """Calculate the similarity of two activity_profiles of the same length.
        When mask is given, only the significant in the mask are taken into account.
        Returns value between 0 .. 1"""

        ref_len = len(reference)
        mask_len = len(mask)
        assert ref_len == len(current), "Activity Profiles must have the same length to be compared."
        matching_bytes = 0
        sig_pos = 0
        for idx in range(ref_len):
            # when a mask is used only consider those positions that are significant in the mask
            if mask_len > 0 and mask[idx] == "1": continue
            if current[idx] != "1" or reference[idx] != "1":
                sig_pos += 1
                if ignore_direction:
                    if current[idx] != "1" and reference[idx] != "1":
                        matching_bytes += 1
                else:
                    if current[idx] == reference[idx]:
                        matching_bytes += 1
        return matching_bytes / sig_pos


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


def print_iter(l):
    for x in l:
        print(x)


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


def middle(lst, size):
    """Return the middle fraction of a sorted list, removing outliers."""
    mid_lst = sorted(list(lst))
    l_mid_lst = len(mid_lst)
    num_el = int(size * l_mid_lst)
    start = (l_mid_lst - num_el) // 2
    end = start + num_el
    mid_lst = mid_lst[start:end]
    return mid_lst


def prof_to_list(act_prof):
    return [int(x) for x in act_prof]


def split_prof(df, id_prop):
    """Split the Activity Profile into individual parameters."""
    props = [id_prop, "Act_Profile"]
    df_props = df[props]
    result = []
    for _, rec in df_props.iterrows():
        ap = rec.pop("Act_Profile")
        for i, v in enumerate(ap):
            rec[ACT_PROF_PARAMETERS[i][7:]] = int(v)
        result.append(rec)
    return pd.DataFrame(result)


def melt(df, id_prop="Compound_Id"):
    """Taken and adapted from the Holoviews measles heatmap example."""
    result = split_prof(df, id_prop)
    result = pd.melt(result, id_vars=id_prop, var_name="Parameter", value_name="Value")
    result = result.reset_index().drop("index", axis=1)
    # result = result.sort_values([id_prop, "Parameter"]).reset_index().drop("index", axis=1)
    result = result[["Parameter", id_prop, "Value"]]
    return result
