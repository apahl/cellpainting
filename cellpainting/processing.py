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
import os.path as op

import pandas as pd
import numpy as np

from . import tools as cpt
from .config import ACT_PROF_PARAMETERS

try:
    from misc_tools import apl_tools as apt
    AP_TOOLS = True
    #: Library version
    VERSION = apt.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))

except ImportError:
    AP_TOOLS = False
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))


FINAL_PARAMETERS = ['Metadata_Plate', 'Metadata_Well', 'plateColumn', 'plateRow',
                    "Compound_Id", 'Batch_Id', "Producer", "KnownAct",
                    'WellType', 'Conc_uM', "Activity", "Act_Profile"]
DROP_FROM_NUMBERS = ['plateColumn', 'plateRow', 'Conc_uM', "Compound_Id"]
DROP_GLOBAL = ["PathName_CellOutlines", "URL_CellOutlines", 'FileName_CellOutlines',
               'ImageNumber', 'Metadata_Site', 'Metadata_Site_1', 'Metadata_Site_2']


class DataSet():
    def __init__(self):
        self.data = pd.DataFrame()
        self.fields = {"plateColumn": "Metadata_Plate",
                       "WellType": "WellType", "ControlWell": "Control", "CompoundWell": "Compound"}


    def xxx_repr_html_(self):
        return self.data._repr_html_()


    def _repr_html_(self):
        parameters = [k for k in FINAL_PARAMETERS if k in self.data]
        print("Shape:     ", self.shape)
        print("Parameters:", parameters)
        return self.data[parameters]._repr_html_()


    def load(self, fn):
        """Read one or multiple result files and concatenate them into one dataset.
        `fn` is a single filename (string) or a list of filenames."""
        self.data = load(fn).data


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
        result = DataSet()
        result.data = well_type_from_position(self.data)
        return result


    def well_from_position(self, well_name="Metadata_Well",
                           row_name="plateRow", col_name="plateColumn"):
        """Assign Metadata_Well from plateRow, plateColumn"""
        result = DataSet()
        result.data = well_from_position(self.data, well_name=well_name,
                                         row_name=row_name, col_name=col_name)
        return result


    def position_from_well(self, well_name="Metadata_Well",
                           row_name="plateRow", col_name="plateColumn"):
        """Generate plateRow and plateColumn from Metatadata_Well"""
        result = DataSet()
        result.data = position_from_well(self.data, well_name=well_name,
                                         row_name=row_name, col_name=col_name)
        return result


    def join_layout_384(self, layout_fn, on="Address"):
        result = DataSet()
        result.data = join_layout_384(self.data, layout_fn, on=on)
        return result


    def join_layout_1536(self, layout_fn, plate, on="Address"):
        """Cell Painting is always run in 384er plates.
        COMAS standard screening plates are format 1536.
        With this function, the 1536-to-384 reformatting file
        with the smiles added by join_smiles_to_layout_1536()
        can be used directly to join the layout to the individual 384er plates."""
        result = DataSet()
        result.data = join_layout_1536(self.data, layout_fn, plate, on=on)
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
        return result, toxic


    def remove_flagged(self):
        """Remove entries with `Pure_Flag == "Fail"`"""
        result = DataSet()
        result.data = remove_flagged(self.data)
        return result


    def remove_outliers(self, times_dev=3.0, group_by=None, method="median"):
        """Returns the filtered dataframe as well as the outliers.
        method can be `median`or `mean` """
        result = DataSet()
        outliers = DataSet()
        result.data, outliers.data = remove_outliers(self.data, times_dev=times_dev,
                                                     group_by=group_by, method=method)
        return result, outliers


    def group_on_well(self, group_by=FINAL_PARAMETERS):
        """Group results on well level."""
        result = DataSet()
        result.data = group_on_well(self.data, group_by=group_by)
        return result


    def poc(self, group_by=None, well_type="WellType", control_name="Control"):
        """Normalize the data set to Percent-Of-Control per group (e.g. per plate)
        based on the median of the controls.
        Parameters:
            group_by (string or None): optional column by which the calculation should be grouped,
            e.g. the column with plate name."""
        result = DataSet()
        result.data = poc(self.data, group_by=group_by)
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
        return result


    def relevant_parameters(self, ctrls_mad_min=0.2, ctrls_mad_max=5.0, times_mad=3.5):
        result = DataSet()
        result.data = relevant_parameters(self.data, ctrls_mad_min=ctrls_mad_min, ctrls_mad_max=ctrls_mad_max,
                                          times_mad=times_mad)
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
        result.data = correlation_filter(self.data, cutoff=cutoff, method=method)
        return result


    def find_similar(self, act_profile, cutoff=0.9):
        """Filter the dataframe for activity profiles similar to the given one.
        `cutoff` gives the similarity threshold, default is 0.9."""
        result = DataSet()
        result.data = find_similar(self.data, act_profile=act_profile, cutoff=cutoff)
        return result


    @property
    def shape(self):
        return self.data.shape


def load(fn):
    """Read one or multiple result files and concatenate them into one dataset.
    `fn` is a single filename (string) or a list of filenames."""
    result = DataSet()
    if isinstance(fn, list):
        result.data = pd.concat((pd.read_csv(f) for f in fn))
    else:
        result.data = pd.read_csv(fn)

    drop = [d for d in DROP_GLOBAL if d in result.data.keys()]
    result.data.drop(drop, axis=1, inplace=True)
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
    result = result.join(layout, on=on)
    result.drop(on, axis=1, inplace=True)
    return result


def join_layout_1536(df, layout_fn, plate, on="Address"):
    """Cell Painting is always run in 384er plates.
    COMAS standard screening plates are format 1536.
    With this function, the 1536-to-384 reformatting file
    with the smiles added by join_smiles_to_layout_1536()
    can be used directly to join the layout to the individual 384er plates."""
    result = df.copy()
    layout = pd.read_csv(layout_fn)
    result[on] = plate[:-1] + result["Metadata_Well"]
    result = result.join(layout, on=on)
    result.drop(on, axis=1, inplace=True)
    return result


def numeric_parameters(df):
    result = df.copy().select_dtypes(include=[np.number])
    drop = [d for d in DROP_FROM_NUMBERS if d in result.keys()]
    result.drop(drop, axis=1, inplace=True)
    return result


def remove_toxic(df, cutoff=0.55):
    """Remove data rows of toxic compounds"""
    median_cell_count_controls = df[df["WellType"] == "Control"]["Count_Cells"].median()
    result = df[df["Count_Cells"] >= median_cell_count_controls * cutoff]
    toxic = df[df["Count_Cells"] < median_cell_count_controls * cutoff]

    return result, toxic


def remove_flagged(df, strict=False, reset_index=True):
    """Remove entries with `Pure_Flag == "Fail"`
    If `strict == True` compound with `Pure_Flag == Warn` are also removed."""
    result = df.copy()
    outliers_list = []
    result = result[result["Pure_Flag"] != "Fail"]
    outl = result[result["Pure_Flag"] == "Fail"]
    outliers_list.append(outl)
    if strict:
        result = result[result["Pure_Flag"] != "Warn"]
        outl = result[result["Pure_Flag"] == "Warn"]
        outliers_list.append(outl)
    outliers = pd.concat(outliers_list)
    if reset_index:
        result = result.reset_index()
        outliers = outliers.reset_index()
    return result, outliers


def _remove_outliers(df, times_mad=3.0):
    include = [k for k in FINAL_PARAMETERS if k in df.keys()]
    result = numeric_parameters(df)
    mask = (result - result.median()).abs() - times_mad * result.mad() <= 0
    result = result[(mask).all(axis=1)]
    # add the non-numeric columns again
    for k in include:
        result[k] = df[k]
    return result


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
        act_parameters = [k for k in df.select_dtypes(include=[np.number]).keys()
                          if k.startswith("Count_") or k.startswith("Mean_")]
        print(len(act_parameters))
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


def relevant_parameters(df, ctrls_mad_min=0.2, ctrls_mad_max=5.0, times_mad=3.5):
    relevant_table = FINAL_PARAMETERS.copy()
    result = df.copy()
    controls = df[df["WellType"] == "Control"].select_dtypes(include=[pd.np.number])
    compounds = df[df["WellType"] == "Compound"].select_dtypes(include=[pd.np.number])

    ds = controls.mad() >= ctrls_mad_min
    ctrl_set = set([p for p in ds.keys() if ds[p]])

    ds = controls.mad() < ctrls_mad_max
    tmp_set = set([p for p in ds.keys() if ds[p]])

    ctrl_set.intersection_update(tmp_set)

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
    print("It took {} iterations to remove all correlated parameters.".format(iteration - 1))
    return df[parameters_uncorr]


def find_similar(df, act_profile, cutoff=0.9):
    """Filter the dataframe for activity profiles similar to the given one.
    `cutoff` gives the similarity threshold, default is 0.9."""
    result = df[df.apply(lambda x: cpt.profile_sim(x["Act_Profile"], act_profile) >= cutoff, axis=1)]
    return result
