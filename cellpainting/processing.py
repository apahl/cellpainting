#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##########
Processing
##########

*Created on Thu Jun  1 14:15 2017 by A. Pahl*

Processing results from the CellPainting Assay in the Jupyter notebook.
This module provides the DataSet class and its methods.
Additional functions in this module act on pandas DataFrames.
convert 170512_B03_s5_w149D3C0B4-85CE-42BE-AB4F-1B460BEECC73.tif -resize 200x200 -negate 170512_B03_s5_w1.png

# invert image:
from PIL import Image
import PIL.ImageOps
image = Image.open('your_image.png')
inverted_image = PIL.ImageOps.invert(image)
inverted_image.save('new_name.png')

# resize an image
size = (200, 200)
im = Image.open(infile)
im.thumbnail(size, Image.ANTIALIAS)
im.save(outfile, "JPEG")
"""

import time
import glob
import os.path as op
from collections import Counter
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np

from IPython.core.display import HTML

from . import tools as cpt
from .config import ACT_PROF_PARAMETERS, ACT_CUTOFF

try:
    from misc_tools import apl_tools, comas_config
    AP_TOOLS = True
    #: Library version
    VERSION = apl_tools.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))
    COMAS = comas_config.COMAS
    ANNOTATIONS = comas_config.ANNOTATIONS
    REFERENCES = comas_config.REFERENCES

except ImportError:
    AP_TOOLS = False
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))
    COMAS = ""
    ANNOTATIONS = ""
    REFERENCES = ""


FINAL_PARAMETERS = ['Metadata_Plate', 'Metadata_Well', 'plateColumn', 'plateRow',
                    "Compound_Id", 'Batch_Id', "Producer", "KnownAct", "Trivial_Name",
                    'WellType', 'Conc_uM', "Activity", "Act_Profile", "Smiles"]
DROP_FROM_NUMBERS = ['plateColumn', 'plateRow', 'Conc_uM', "Compound_Id"]
DROP_GLOBAL = ["PathName_CellOutlines", "URL_CellOutlines", 'FileName_CellOutlines',
               'ImageNumber', 'Metadata_Site', 'Metadata_Site_1', 'Metadata_Site_2']


class DataSet():
    def __init__(self, log=True):
        self.data = pd.DataFrame()
        self.fields = {"plateColumn": "Metadata_Plate",
                       "WellType": "WellType", "ControlWell": "Control", "CompoundWell": "Compound"}
        self.log = log


    def __getitem__(self, item):
        result = DataSet()
        result.data = self.data[item]
        return result


    def __getattr__(self, name):
        """Try to call undefined methods on the underlying pandas DataFrame."""
        def method(*args, **kwargs):
            return getattr(self.data, name)(*args, **kwargs)
        return method


    def show(self):
        parameters = [k for k in FINAL_PARAMETERS if k in self.data]
        print("Shape:     ", self.shape)
        print("Parameters:", parameters)
        return HTML(self.data[parameters]._repr_html_())


    def head(self, n=5):
        parameters = [k for k in FINAL_PARAMETERS if k in self.data]
        print("Shape:     ", self.shape)
        return self.data[parameters].head(n)


    def drop_cols(self, cols, inplace=False):
        if inplace:
            self.data.drop(cols, axis=1, inplace=True)
            self.print_log("drop cols (inplace)")
        else:
            result = DataSet()
            result.data = self.data.drop(cols, axis=1, inplace=False)
            result.print_log("drop cols")
            return result


    def keep_cols(self, cols, inplace=False):
        if inplace:
            self.data = self.data[cols]
            self.print_log("keep cols (inplace)")
        else:
            result = DataSet()
            result.data = self.data[cols]
            result.print_log("keep cols")
            return result


    def print_log(self, component, add_info=""):
        if self.log:
            component = component + ":"
            if len(add_info) > 0:
                add_info = "    ({})".format(add_info)
            print("{:22s} ({:4d} | {:4d}){}".format(component, self.shape[0], self.shape[1], add_info))


    def load(self, fn, sep="\t"):
        """Read one or multiple result files and concatenate them into one dataset.
        `fn` is a single filename (string) or a list of filenames."""
        self.data = load(fn, sep=sep).data
        self.print_log("load data")


    def write_csv(self, fn, parameters=None, sep="\t"):
        result = self.data.copy()
        if isinstance(parameters, list):
            result = result[parameters]
        result.to_csv(fn, sep=sep, index=False)


    def describe(self, times_mad=3.0):
        df = numeric_parameters(self.data)
        stats = pd.DataFrame()
        stats["Min"] = df.min()
        stats["Max"] = df.max()
        stats["Median"] = df.median()
        stats["MAD"] = df.mad()
        stats["Outliers"] = df[(((df - df.median()).abs() - times_mad * df.mad()) > 0)].count()
        print(self.shape)
        return stats


    def well_type_from_position(self):
        """Assign the WellType from the position on the plate.
        Controls are in column 11 and 12"""
        result = DataSet(log=self.log)
        result.data = well_type_from_position(self.data)
        return result


    def well_from_position(self, well_name="Metadata_Well",
                           row_name="plateRow", col_name="plateColumn"):
        """Assign Metadata_Well from plateRow, plateColumn"""
        result = DataSet(log=self.log)
        result.data = well_from_position(self.data, well_name=well_name,
                                         row_name=row_name, col_name=col_name)
        return result


    def position_from_well(self, well_name="Metadata_Well",
                           row_name="plateRow", col_name="plateColumn"):
        """Generate plateRow and plateColumn from Metatadata_Well"""
        result = DataSet(log=self.log)
        result.data = position_from_well(self.data, well_name=well_name,
                                         row_name=row_name, col_name=col_name)
        return result


    def join_layout_384(self, layout_fn, on="Address_384"):
        result = DataSet(log=self.log)
        result.data = join_layout_384(self.data, layout_fn, on=on)
        result.print_log("join layout 384")
        return result


    def join_layout_1536(self, layout_fn, plate, on="Address_384"):
        """Cell Painting is always run in 384er plates.
        COMAS standard screening plates are format 1536.
        With this function, the 1536-to-384 reformatting file
        with the smiles added by join_smiles_to_layout_1536()
        can be used directly to join the layout to the individual 384er plates."""
        result = DataSet(log=self.log)
        result.data = join_layout_1536(self.data, layout_fn, plate, on=on)
        result.print_log("join layout 1536")
        return result


    def numeric_parameters(self):
        result = DataSet()
        result.data = numeric_parameters(self.data)
        return result


    def remove_toxic(self, cutoff=0.55):
        """Remove data rows of toxic compounds"""
        result = DataSet()
        toxic = DataSet()
        result.data, toxic.data = remove_toxic(self.data, cutoff=cutoff)
        result.print_log("remove toxic", "{:3d} removed".format(toxic.shape[0]))
        return result, toxic


    def remove_flagged(self, strict=False, reset_index=True):
        """Remove entries with `Pure_Flag == "Fail"`"""
        result = DataSet()
        flagged = DataSet()
        result.data, flagged.data = remove_flagged(self.data)
        result.print_log("remove flagged", "{:3d} removed".format(flagged.shape[0]))
        return result, flagged


    def remove_outliers(self, times_dev=3.0, group_by=None, method="median"):
        """Returns the filtered dataframe as well as the outliers.
        method can be `median`or `mean` """
        result = DataSet()
        outliers = DataSet()
        result.data, outliers.data = remove_outliers(self.data, times_dev=times_dev,
                                                     group_by=group_by, method=method)
        result.print_log("remove outliers", "{:3d} removed".format(outliers.shape[0]))
        return result, outliers


    def remove_skipped_echo_direct_transfer(self, fn):
        """Remove wells that were reported as skipped in the Echo protocol (xml).
        This functions works with Echo direct transfer protocols.
        Function supports using wildcards in the filename, the first file will be used.
        Returns a new dataframe without the skipped wells."""
        result = DataSet()
        result.data, skipped = remove_skipped_echo_direct_transfer(self.data, fn=fn)
        skipped_str = "(" + ", ".join(skipped) + ")"
        result.print_log("remove skipped", "{:3d} skipped {}".format(self.shape[0] - result.shape[0],
                                                                     skipped_str))
        return result


    def drop_dups(self, cpd_id="Compound_Id"):
        """Drop duplicate Compound_Ids"""
        result = DataSet()
        result.data = self.data.drop_duplicates(cpd_id)
        result.print_log("drop dups")
        return result


    def group_on_well(self, group_by=FINAL_PARAMETERS):
        """Group results on well level."""
        result = DataSet()
        result.data = group_on_well(self.data, group_by=group_by)
        result.print_log("group on well")
        return result


    def join_smiles(self, df_smiles=None):
        """Join Smiles from Compound_Id."""
        result = DataSet()
        result.data = join_smiles(self.data, df_smiles=df_smiles)
        result.print_log("join smiles")
        return result


    def join_annotations(self):
        """Join Annotations from Compound_Id."""
        result = DataSet()
        result.data = join_annotations(self.data)
        result.print_log("join annotations")
        return result


    def poc(self, group_by=None, well_type="WellType", control_name="Control"):
        """Normalize the data set to Percent-Of-Control per group (e.g. per plate)
        based on the median of the controls.
        Parameters:
            group_by (string or None): optional column by which the calculation should be grouped,
            e.g. the column with plate name."""
        result = DataSet()
        result.data = poc(self.data, group_by=group_by)
        self.print_log("POC")
        return result


    def activity_profile(self, mad_mult=3.5, parameters=ACT_PROF_PARAMETERS, only_final=True):
        """Generates the `Act_Profile` column.
        The byte is set when the parameter's value is greater (or smaller)
        than parameter_ctrl.median() + (or -) `mad_mult`* parameter.mad()

        If a list of parameters is given, then the activity profile will be calculated for these parameters.

        If `only_final` == `True`, then only the parameters listed in `FINAL_PARAMETERS` are kept in the output_table.

        Returns a new Pandas DataFrame."""
        result = DataSet()
        result.data = activity_profile(self.data, mad_mult=mad_mult, parameters=parameters,
                                       only_final=only_final)
        result.print_log("activity profile")
        return result


    def relevant_parameters(self, ctrls_mad_rel_min=0.01, ctrls_mad_rel_max=0.10, times_mad=3.5):
        result = DataSet()
        result.data = relevant_parameters(self.data, ctrls_mad_rel_min=ctrls_mad_rel_min,
                                          ctrls_mad_rel_max=ctrls_mad_rel_max, times_mad=times_mad)
        result.print_log("relevant parameters")
        return result


    def correlation_filter(self, cutoff=0.9, method="pearson"):
        """The correlation removes all highly correlated columns from the dataframe.
        The function was implemented according to the description of the corresponding
        KNIME component.

        Parameters:
            cutoff (float): correlation cutoff
            method (string): "pearson", "kendall", "spearman" (very slow)

        Returns a new DataFrame with only the non-correlated columns"""
        result = DataSet()
        result.data, iterations = correlation_filter(self.data, cutoff=cutoff, method=method)
        result.print_log("correlation filter", "{:3d} iterations".format(iterations))
        return result


    def add_act_profile_for_control(self, parameters=ACT_PROF_PARAMETERS):
        # Compound_Id DMSO: 245754
        control = {"Compound_Id": 245754, "Trivial_Name": "Control", "Activity": 0,
                   "Act_Profile": "".join(["1"] * len(parameters))}
        ck = control.keys()
        for k in ck:
            if k not in self.data.keys():
                control.pop(k)
        tmp = pd.DataFrame(control)
        result = DataSet()
        result.data = pd.concat(self.data, tmp)
        return result


    def find_similar(self, act_profile, cutoff=0.9, max_num=5):
        """Filter the dataframe for activity profiles similar to the given one.
        `cutoff` gives the similarity threshold, default is 0.9."""
        result = DataSet()
        result.data = find_similar(self.data, act_profile=act_profile, cutoff=cutoff, max_num=max_num)
        result.print_log("find similar")
        return result


    def find_similar_in_refs(self, act_profile, df_refs=None, cutoff=0.9, sep="\t"):
        """Find and add references with similar activity profiles to the dataframe
        `cutoff` gives the similarity threshold, default is 0.9."""
        result = DataSet()
        result.data = find_similar_in_refs(self.data, act_profile=act_profile,
                                           cutoff=cutoff, sep=sep)
        # result.print_log("find similar")
        return result


    def count_active_parameters_occurrences(self, parameters=ACT_PROF_PARAMETERS):
        """Counts the number of times each parameter has been active in the dataset."""
        return count_active_parameters_occurrences(self.data, parameters=ACT_PROF_PARAMETERS)


    @property
    def shape(self):
        return self.data.shape


    @property
    def metadata(self):
        """Returns a list of the those parameters in the DataFrame that are NOT CellProfiler measurements."""
        return metadata(self.data)


    @property
    def measurements(self):
        """Returns a list of the CellProfiler parameters that are in the DataFrame."""
        return measurements(self.data)


def load(fn, sep="\t"):
    """Read one or multiple result files and concatenate them into one dataset.
    `fn` is a single filename (string) or a list of filenames."""
    result = DataSet()
    if isinstance(fn, list):
        result.data = pd.concat((pd.read_csv(f, sep=sep) for f in fn))
    else:
        result.data = pd.read_csv(fn, sep=sep)

    drop = [d for d in DROP_GLOBAL if d in result.data.keys()]
    result.data.drop(drop, axis=1, inplace=True)
    result.print_log("load dataset")
    return result


def read_smiles_file(fn, props=['Compound_Id', "Producer", "Smiles", "Pure_Flag"]):
    """Read in the file with the Compound_Ids and the Simles.
    Return a DataFrame for fast access."""
    result = pd.read_csv(fn, sep="\t")
    result = result[props]
    return result


def well_type_from_position(df):
    """Assign the WellType from the position on the plate.
    Controls are in column 11 and 12"""
    result = df.copy()
    result["WellType"] = "Compound"
    result["WellType"][(result["plateColumn"] == 11) | (result["plateColumn"] == 12)] = "Control"
    return result


def well_from_position(df, well_name="Metadata_Well",
                       row_name="plateRow", col_name="plateColumn"):
    """Assign Metadata_Well from plateRow, plateColumn"""
    def _well_from_position_series(s):
        return cpt.well_from_position_single(s[0], s[1])

    result = df.copy()
    result[well_name] = result[[row_name, col_name]].apply(_well_from_position_series, axis=1)
    return result


def position_from_well(df, well_name="Metadata_Well",
                       row_name="plateRow", col_name="plateColumn"):
    """Generate plateRow and plateColumn from Metatadata_Well"""
    def _position_from_well_series(well):
        return(pd.Series(cpt.position_from_well_single(well)))

    result = df.copy()
    result[[row_name, col_name]] = result[well_name].apply(_position_from_well_series)
    return result


def join_layout_384(df, layout_fn, on="Address"):
    result = df.copy()
    result[on] = result["Metadata_Well"]
    layout = pd.read_csv(layout_fn)
    result = result.merge(layout, on=on)
    result.drop(on, axis=1, inplace=True)
    result = result.apply(pd.to_numeric, errors='ignore')
    return result


def join_layout_1536(df, layout_fn, plate, on="Address_384", sep="\t"):
    """Cell Painting is always run in 384er plates.
    COMAS standard screening plates are format 1536.
    With this function, the 1536-to-384 reformatting file
    with the smiles added by join_smiles_to_layout_1536()
    can be used directly to join the layout to the individual 384er plates."""
    result = df.copy()
    layout = pd.read_csv(layout_fn, sep=sep)
    layout[on] = layout["Plate_name_384"].str[-1:] + layout[on]
    layout["Batch_Id"] = layout["Container_ID_1536"].str[:9]
    layout["Compound_Id"] = layout["Container_ID_1536"].str[:6]
    drop = ["Container_ID_1536", "Plate_name_384", "Plate_name_1536", "Address_1536", "Index"]
    layout.drop(drop, axis=1, inplace=True)
    result[on] = plate[-1:] + result["Metadata_Well"]
    result = result.merge(layout, on=on, how="inner")
    result.drop(on, axis=1, inplace=True)
    result = result.apply(pd.to_numeric, errors='ignore')
    return result


def join_smiles(df, df_smiles=None):
    """Join Smiles from Compound_Id."""
    keep = ['Compound_Id', "Producer", "Smiles", "Pure_Flag"]
    # df["Compound_Id"] = df["Compound_Id"].astype(int)
    if df_smiles is None:
        df_smiles = pd.read_csv(COMAS, sep="\t")
        df_smiles = df_smiles[keep]
    elif isinstance(df_smiles, str):
        df_smiles = pd.read_csv(df_smiles, sep="\t")
    result = df.merge(df_smiles, on="Compound_Id", how="inner")
    result = result.apply(pd.to_numeric, errors='ignore')
    return result


def join_annotations(df):
    """Join Annotations from Compound_Id."""
    annotations = pd.read_csv(ANNOTATIONS, sep="\t")
    result = df.merge(annotations, on="Compound_Id", how="left")
    result = result.replace(np.nan, "", regex=True)
    return result


def metadata(df):
    """Returns a list of the those parameters in the DataFrame that are NOT CellProfiler measurements."""
    parameters = [k for k in df.keys()
                  if not (k.startswith("Count_") or k.startswith("Mean_"))]
    return parameters


def measurements(df):
    """Returns a list of the CellProfiler parameters that are in the DataFrame."""
    parameters = [k for k in df.select_dtypes(include=[np.number]).keys()
                  if k.startswith("Count_") or k.startswith("Mean_")]
    return parameters


def numeric_parameters(df):
    result = df.copy()[measurements(df)]
    return result


def remove_toxic(df, cutoff=0.55):
    """Remove data rows of toxic compounds"""
    median_cell_count_controls = df[df["WellType"] == "Control"]["Count_Cells"].median()
    result = df[df["Count_Cells"] >= median_cell_count_controls * cutoff]
    toxic = df[df["Count_Cells"] < median_cell_count_controls * cutoff]
    return result, toxic


def remove_skipped_echo_direct_transfer(df, fn):
    """Remove wells that were reported as skipped in the Echo protocol (xml).
    This functions works with Echo direct transfer protocols.
    Function supports using wildcards in the filename, the first file will be used.
    Returns a new dataframe without the skipped wells."""
    assert fn.endswith(".xml"), "Echo file expected in XML format."
    skipped_wells = []
    echo_fn = glob.glob(fn)[0]  # use the first glob match
    echo_print = ET.parse(echo_fn).getroot()
    skipped = echo_print.find("skippedwells")
    for well in skipped.findall("w"):
        skipped_wells.append(cpt.format_well(well.get("dn")))
    # print("Skipped wells (will be removed):", skipped_wells)
    # remove the rows with the skipped wells
    #   i.e. keep the rows where Metadata_Well is not in the list skipped_wells
    result = df[~df["Metadata_Well"].isin(skipped_wells)]
    return result, skipped_wells


def remove_flagged(df, strict=False, reset_index=True):
    """Remove entries with `Pure_Flag == "Fail"`
    If `strict == True` compound with `Pure_Flag == Warn` are also removed."""
    result = df.copy()
    outliers_list = []
    outl = result[result["Pure_Flag"] == "Fail"]
    result = result[result["Pure_Flag"] != "Fail"]
    print(outl.shape)
    outliers_list.append(outl)
    if strict:
        outl = result[result["Pure_Flag"] == "Warn"]
        result = result[result["Pure_Flag"] != "Warn"]
        outliers_list.append(outl)
    outliers = pd.concat(outliers_list)
    if reset_index:
        result = result.reset_index()
        outliers = outliers.reset_index()
        result.drop("index", axis=1, inplace=True)
        outliers.drop("index", axis=1, inplace=True)
    return result, outliers


def remove_outliers(df, times_dev=3.0, group_by=None, method="median", reset_index=True):
    """Returns the filtered dataframe as well as the outliers.
    method can be `median`or `mean` """
    include = [k for k in FINAL_PARAMETERS if k in df.keys()]
    input = df.copy()
    # input = numeric_parameters(df)
    if group_by is None:
        group_by = "temp_group"
        input[group_by] = "data"

    gdata_list = []
    outliers_list = []
    for group in df[group_by].unique():
        gdata = input[input[group_by] == group]
        gdata = numeric_parameters(gdata)
        if method == "median":
            mask = (gdata - gdata.median()).abs() - times_dev * gdata.mad() <= 0
        elif method == "mean":
            mask = (gdata - gdata.mean()).abs() - times_dev * gdata.std() <= 0
        else:
            raise ValueError("Unknown method {}.".format(method))
        good_data = gdata[(mask).all(axis=1)]
        outl_data = gdata[(~(mask).all(axis=1))]  # outliers
        # print(group, ": ", good_data.shape, outl_data.shape)
        gdata_list.append(good_data)
        outliers_list.append(outl_data)
    result = pd.concat(gdata_list)
    outliers = pd.concat(outliers_list)

    if group_by == "temp_group":  # remove the grouping temp col again
        result.drop(group_by, axis=1, inplace=True)
        outliers.drop(group_by, axis=1, inplace=True)
    for k in include:
        result[k] = df[k]
        outliers[k] = df[k]
    if reset_index:
        result = result.reset_index()
        outliers = outliers.reset_index()
    return result, outliers


def group_on_well(df, group_by=FINAL_PARAMETERS):
    """Group results on well level."""
    group_by = list(set(group_by).intersection(set(df.keys())))
    result = df.groupby(by=group_by).median().reset_index()
    return result


def poc(df, group_by=None):
    result = df.copy()
    if group_by is None:  # create a temp grouping column
        group_by = "temp_group"
        result[group_by] = "data"

    plates = set(result[group_by])
    for plate in plates:
        print("Normalizing {}...   ".format(plate), end="")
        controls = result[(result[group_by] == plate) & (result["WellType"] == "Control")].select_dtypes(include=[np.number])
        median = controls.median()
        for col in controls.keys():
            if col in FINAL_PARAMETERS: continue
            result[col] = 100 * result[col] / median[col]
        print("done.")

    if group_by == "temp_group":  # remove the grouping temp col again
        result.drop(group_by, axis=1, inplace=True)
    return result


def activity_profile(df, mad_mult=3.5, parameters=ACT_PROF_PARAMETERS, only_final=True):
    """Generates the `Act_Profile` column.
    The byte is set when the parameter's value is greater (or smaller)
    than parameter_ctrl.median() + (or -) `mad_mult`* parameter.mad()

    If a list of parameters is given, then the activity profile will be calculated for these parameters.

    If `only_final` == `True`, then only the parameters listed in `FINAL_PARAMETERS`
    are kept in the output_table.

    Returns a new Pandas DataFrame."""
    result = df.copy()

    if parameters is None:  # choose all numeric parameters
        act_parameters = measurements(df)
    else:
        act_parameters = parameters.copy()
    # sort parameters alphabetically
    act_parameters.sort()
    controls = df[act_parameters][df["WellType"] == "Control"]

    for key in act_parameters:
        median = controls[key].median()
        times_mad = mad_mult * controls[key].mad()
        lower_bound = median - times_mad
        upper_bound = median + times_mad
        result.loc[df[key].between(lower_bound, upper_bound, inclusive=True), [key]] = 1
        result.loc[df[key] < lower_bound, [key]] = 0
        result.loc[df[key] > upper_bound, [key]] = 2

    result[act_parameters] = result[act_parameters].astype(int)
    result["Activity"] = (result[act_parameters] != 1).sum(axis=1)
    result["Act_Profile"] = result[act_parameters].astype(str).apply(lambda x: "".join(x), axis=1)

    if only_final:
        for k in result.keys():
            if k not in FINAL_PARAMETERS:
                result.drop(k, axis=1, inplace=True)
    return result


def relevant_parameters(df, ctrls_mad_rel_min=0.01, ctrls_mad_rel_max=0.10, times_mad=3.5):
    """...mad_rel...: MAD relative to the median value"""
    relevant_table = FINAL_PARAMETERS.copy()
    result = df.copy()
    controls = df[df["WellType"] == "Control"].select_dtypes(include=[pd.np.number])
    compounds = df[df["WellType"] == "Compound"].select_dtypes(include=[pd.np.number])

    ds = controls.mad() / controls.median() >= ctrls_mad_rel_min
    ctrl_set = set([p for p in ds.keys() if ds[p]])
    print(len(ctrl_set))

    ds = controls.mad() / controls.median() < ctrls_mad_rel_max
    tmp_set = set([p for p in ds.keys() if ds[p]])
    print(len(tmp_set))

    ctrl_set.intersection_update(tmp_set)
    print(len(ctrl_set))

    controls = controls[list(ctrl_set)]
    compounds = compounds[list(ctrl_set)]

    ds = compounds.max() - controls.median() - times_mad * controls.mad() > 0
    cpd_max_set = set([p for p in ds.keys() if ds[p]])
    ds = controls.median() - compounds.min() - times_mad * controls.mad() > 0
    cpd_min_set = set([p for p in ds.keys() if ds[p]])

    cpd_set = cpd_max_set.union(cpd_min_set)

    for p in cpd_set:
        relevant_table.append(p)

    result_keys = list(result.keys())
    for key in result_keys:
        if key not in relevant_table:
            result.drop(key, axis=1, inplace=True)
    return result


def correlation_filter(df, cutoff=0.9, method="pearson"):
    assert method in ["pearson", "kendall", "spearman"], 'method has to be one of ["pearson", "kendall", "spearman"]'

    df_copy = df.copy().select_dtypes(include=[np.number])

    # init the list of the uncorrelated parameters, incl. some string param.
    parameters_uncorr = [p for p in FINAL_PARAMETERS if p in df]

    iteration = 0
    while True:
        cm = df_copy.corr(method=method)

        ds = cm[cm > cutoff].count().sort_values(ascending=False)
        if ds[0] == 1: break  # no more correlations
        iteration += 1

        # from all columns with the same number of correlated columns,
        # find the column with the highest value range
        # and keep that preferably
        num_correlated = ds[0]  # number of correlated columns
        rnge = 0.0
        rnge_key = ""
        for i in range(len(ds)):
            if ds[i] < num_correlated: break  # only compare columns with the same highest correlation
            k = ds.keys()[i]
            r = df_copy[k].max() - df_copy[k].min()
            if r > rnge:
                rnge = r
                rnge_key = k

        # print(rnge_key, "  (", num_correlated, ")")
        keep_it = rnge_key
        # print("keep_it:", keep_it)
        parameters_uncorr.append(keep_it)

        parameters_to_remove = list(cm[keep_it][cm[keep_it] > cutoff].keys())
        # The uncorrelated parameter `keep_it` is also in this list and has to be removed from it:
        # parameters_to_remove.remove(keep_it)
        # print("parameters_to_remove:", parameters_to_remove)

        # remove the correlated parameters from both axes:
        df_copy.drop(parameters_to_remove, axis=1, inplace=True)

    parameters_uncorr.extend(df_copy.keys())
    parameters_uncorr = list(set(parameters_uncorr))
    # print("It took {} iterations to remove all correlated parameters.".format(iteration - 1))
    return df[parameters_uncorr], iteration


def find_similar(df, act_profile, cutoff=0.9, max_num=5):
    """Filter the dataframe for activity profiles similar to the given one.
    `cutoff` gives the similarity threshold, default is 0.9."""
    result = df.copy()
    result["Similarity"] = result["Act_Profile"].apply(lambda x: cpt.profile_sim(x["Act_Profile"],
                                                                                 act_profile))
    result = result[result["Similarity"] >= cutoff]
    result.drop("Act_Profile", axis=1, inplace=True)
    result = result.sort_values("Similarity", ascending=False).head(max_num)
    return result


def find_similar_in_refs(df, act_profile, df_refs=None, cutoff=0.9, sep="\t"):
    """Find and add references with similar activity profiles to the dataframe
    `cutoff` gives the similarity threshold, default is 0.9."""
    if df_refs is None:
        df_refs = load(REFERENCES, sep=sep)
    elif isinstance(df_refs, str):
        df_refs = load(df_refs, sep=sep)
    df_refs.log = False
    result = {"Compound_Id": [], "Activity": [], "Similar_Ref": [], "Similarity": [],
              "Trivial_Name": [], "Annotation": [], "Act_Flag": []}
    for _, rec in df.iterrows():
        if rec["Activity"] < ACT_CUTOFF:
            result["Compound_Id"].append(rec["Compound_Id"])
            result["Activity"].append(rec["Activity"])
            result["Act_Flag"].append("InAct")
            for k in ["Similar_Ref", "Similarity", "Annotation", "Trivial_Name"]:
                result[k].append("")
        else:
            act_prof = rec["Activity_Profile"]
            sim_refs = df_refs.find_similar(act_prof)
            for _, ref in sim_refs.iterrows():
                result["Compound_Id"].append(rec["Compound_Id"])
                result["Activity"].append(rec["Activity"])
                result["Act_Flag"].append("Active")
                result["Similar_Ref"].append(ref["Compound_Id"])
                result["Similarity"].append(ref["Similarity"])
                result["Trivial_Name"].append(ref["Trivial_Name"])
                result["Annotation"].append(ref["Annotation"])
    return result


def count_active_parameters_occurrences(df, parameters=ACT_PROF_PARAMETERS):
    """Counts the number of times each parameter has been active in the dataset."""
    ctr_int = Counter()
    ctr_str = Counter()
    for _, rec in df.iterrows():
        for idx, b in enumerate(rec["Act_Profile"]):
            if b != "1":
                ctr_int[idx] += 1
    for idx in ctr_int:
        ctr_str[parameters[idx]] = ctr_int[idx]
    return ctr_str
