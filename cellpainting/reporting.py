#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#########
Reporting
#########

*Created on Thu Jun  8 14:40 2017 by A. Pahl*

Tools for creating HTML Reports."""

import time
import base64
import os.path as op
from string import Template
from io import BytesIO as IO

# import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

from PIL import Image, ImageChops

from . import tools as cpt
from . import report_templ as cprt
from . import processing as cpp
from .config import ACT_CUTOFF_PERC

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

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

try:
    # Try to import Avalon so it can be used for generation of 2d coordinates.
    from rdkit.Avalon import pyAvalonTools as pyAv
    USE_AVALON_2D = True
except ImportError:
    print("  * Avalon not available. Using RDKit for 2d coordinate generation.")
    USE_AVALON_2D = False


def check_2d_coords(mol, force=False):
    """Check if a mol has 2D coordinates and if not, calculate them."""
    if not force:
        try:
            mol.GetConformer()
        except ValueError:
            force = True  # no 2D coords... calculate them

    if force:
        if USE_AVALON_2D:
            pyAv.Generate2DCoords(mol)
        else:
            mol.Compute2DCoords()


def mol_from_smiles(smi, calc_2d=True):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles("C")
    else:
        if calc_2d:
            check_2d_coords(mol)
    return mol


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    return None  # no contents


def get_value(str_val):
    if not str_val:
        return ""
    try:
        val = float(str_val)
        if "." not in str_val:
            val = int(val)
    except ValueError:
        val = str_val
    return val


def isnumber(x):
    """Returns True, if x is a number (i.e. can be converted to float)."""
    try:
        float(x)
        return True
    except:
        return False


def convert_bool(dict, dkey, true="Yes", false="No", default="n.d."):
    if dkey in dict:
        if dict[dkey]:
            dict[dkey] = true
        else:
            dict[dkey] = false
    else:
        dict[dkey] = default


def b64_img(mol, size=300):
    img_file = IO()
    img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img.save(img_file, format='PNG')

    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()

    return b64


def mol_img_tag(mol):
    img_tag = '<img src="data:image/png;base64,{}" alt="Mol"/>'.format(b64_img(mol))
    return img_tag


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


def write_page(page, title="Report", fn="index.html"):
    t = Template(cprt.HTML_INTRO + page + cprt.HTML_EXTRO)
    result = t.substitute(title=title)
    write(result, fn=fn)


def overview_report(df, df_refs=None, cutoff=0.6, detailed_cpds=None):
    """detailed_cpds = {Compound_Id: [{"Sim_Ref_Id": Sim_Ref_Id with highest Similarity, "Trivial_Name": ...,
    "Similarity": ..., "Known_Act": ..., "Smiles": ...},
    {{"Sim_Ref_Id": Sim_Ref_Id with next highest Similarity, "Tr...}]}"""
    if isinstance(df, cpp.DataSet):
        df = df.data
    df_refs = cpt.check_df(df_refs, REFERENCES)
    report = [cprt.TABLE_INTRO, cprt.OVERVIEW_TABLE_HEADER]
    row_templ = Template(cprt.OVERVIEW_TABLE_ROW)
    for idx, rec in df.iterrows():
        mol = mol_from_smiles(rec.get("Smiles", "C"))
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx + 1
        convert_bool(rec, "Toxic")

        if rec["Activity"] < ACT_CUTOFF_PERC:
            rec["Act_Flag"] = "inactive"
            rec["Num_Sim_Ref"] = ""
            rec["Link"] = ""
        else:
            rec["Act_Flag"] = "active"
            act_profile = rec["Act_Profile"]
            sim_refs = cpp.find_similar(df_refs, act_profile, cutoff=cutoff)
            if sim_refs.shape[0] == 0:
                rec["Num_Sim_Ref"] = "None"
                rec["Link"] = ""
            else:
                rec["Num_Sim_Ref"] = str(sim_refs.shape[0])
                rec["Link"] = '<a href="details/{}">Detailed<br>Report</a>'.format(rec["Compound_Id"])
                if detailed_cpds is not None:
                    detailed_cpds[rec["Compound_Id"]] = []
                    ref_dict = {}
                    for _, ref in sim_refs.iterrows():
                        ref_dict["Sim_Ref_Id"] = ref["Compound_Id"]
                        ref_dict["Smiles"] = ref["Smiles"]
                        ref_dict["Trivial_Name"] = ref["Trivial_Name"]
                        ref_dict["Similarity"] = ref["Similarity"]
                        ref_dict["Known_Act"] = ref["Known_Act"]
        report.append(row_templ.substitute(rec))
    report.append(cprt.TABLE_EXTRO)
    return "\n".join(report)


def full_report(df, df_refs=None, dirname="report", plate=None, cutoff=0.6):
    cpt.create_dirs(op.join(dirname, "details"))
    if isinstance(df, cpp.DataSet):
        df = df.data
    df_refs = cpt.check_df(df_refs, REFERENCES)
    detailed_cpds = {}
    print("* creating overview...")
    header = "<h2>Cell Painting Overview Report</h2>\n"
    title = "Overview"
    if plate is not None:
        title = plate
        header += "<h3>Plate {}</h3>\n".format(plate)
    header += "<h3>({})</h3>\n".format(time.strftime("%d-%m-%Y %H:%M", time.localtime()))
    overview = header + overview_report(df, df_refs, cutoff=cutoff, detailed_cpds=detailed_cpds)
    write_page(overview, title=title, fn="report/index.html")
