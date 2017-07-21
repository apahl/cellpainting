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
from .config import ACT_PROF_PARAMETERS, ACT_CUTOFF_PERC

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from IPython.core.display import HTML

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
    """detailed_cpds = {Compound_Id: [{"Compound_Id": Sim_Ref_Id with highest Similarity, "Trivial_Name": ...,
    "Similarity": ..., "Known_Act": ..., "Smiles": ...},
    {{"Compound_Id": Sim_Ref_Id with next highest Similarity, "Tr...}]}"""
    if isinstance(df, cpp.DataSet):
        df = df.data
    df_refs = cpt.check_df(df_refs, REFERENCES)
    report = [cprt.TABLE_INTRO, cprt.OVERVIEW_TABLE_HEADER]
    row_templ = Template(cprt.OVERVIEW_TABLE_ROW)
    idx = 0
    for _, rec in df.iterrows():
        idx += 1
        mol = mol_from_smiles(rec.get("Smiles", "C"))
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx
        convert_bool(rec, "Toxic")

        if rec["Activity"] < ACT_CUTOFF_PERC:
            rec["Act_Flag"] = "inactive"
            # rec["Num_Sim_Ref"] = ""
            rec["Max_Sim"] = ""
            rec["Link"] = ""
        else:
            rec["Act_Flag"] = "active"
            act_profile = rec["Act_Profile"]
            sim_refs = cpp.find_similar(df_refs, act_profile, cutoff=cutoff, max_num=5)
            sim_refs = sim_refs.fillna("&mdash;")
            if sim_refs.shape[0] == 0:
                # rec["Num_Sim_Ref"] = "None"
                rec["Max_Sim"] = "< {}".format(cutoff * 100)
            else:
                # rec["Num_Sim_Ref"] = str(sim_refs.shape[0])
                rec["Max_Sim"] = str(sim_refs["Similarity"].max() * 100)
            rec["Link"] = '<a href="details/{}.html">Detailed<br>Report</a>'.format(rec["Compound_Id"])
            if detailed_cpds is not None:
                detailed_cpds[rec["Compound_Id"]] = []
                for _, ref in sim_refs.iterrows():
                    ref_dict = {}
                    ref_dict["Compound_Id"] = ref["Compound_Id"]
                    ref_dict["Smiles"] = ref["Smiles"]
                    ref_dict["Trivial_Name"] = ref["Trivial_Name"]
                    ref_dict["Similarity"] = ref["Similarity"]
                    ref_dict["Known_Act"] = ref["Known_Act"]
                    detailed_cpds[rec["Compound_Id"]].append(ref_dict)
        report.append(row_templ.substitute(rec))
    report.append(cprt.TABLE_EXTRO)
    return "\n".join(report)


def sim_ref_table(sim_refs):
    table = []
    templ = Template(cprt.REF_TABLE_ROW)
    for idx, sim_ref in enumerate(sim_refs, 1):
        rec = sim_ref.copy()
        mol = mol_from_smiles(rec.get("Smiles", "C"))
        rec["Similarity"] *= 100
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx
        row = templ.substitute(rec)
        table.append(row)
    return "\n".join(table)


def changed_parameters_table(act_prof, val, parameters=ACT_PROF_PARAMETERS):
    changed = cpt.parameters_from_act_profile_by_val(act_prof, val, parameters=parameters)
    table = []
    templ = Template(cprt.PARM_TABLE_ROW)
    for idx, p in enumerate(changed, 1):
        rec = {
            "idx": idx,
            "Parameter": p[5:]
        }
        row = templ.substitute(rec)
        table.append(row)
    return "\n".join(table)


def detailed_report(rec, sim_refs):
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
    act_prof = rec["Act_Profile"]
    inc_parm = changed_parameters_table(act_prof, "2")
    dec_parm = changed_parameters_table(act_prof, "0")
    mol = mol_from_smiles(rec.get("Smiles", "C"))

    templ_dict = {
        "Compound_Id": rec["Compound_Id"],
        "Date": date,
        "Producer": rec["Producer"],
        "Activity": rec["Activity"],
        "mol_img": mol_img_tag(mol),
        "Inc_Parm_Table": inc_parm,
        "Dec_Parm_Table": dec_parm
    }
    if len(sim_refs) > 0:
        ref_tbl = sim_ref_table(sim_refs)
        templ_dict["Ref_Table"] = ref_tbl
        t = Template(cprt.DETAILS_TEMPL)
    else:
        t = Template(cprt.DETAILS_TEMPL_NO_SIM_REFS)

    report = t.substitute(templ_dict)
    return report


def full_report(df, df_refs=None, dirname="report", plate=None, cutoff=0.6):
    overview_fn = op.join(dirname, "index.html")
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
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
    header += "<p>({})</p>\n".format(date)
    overview = header + overview_report(df, df_refs, cutoff=cutoff, detailed_cpds=detailed_cpds)
    write_page(overview, title=title, fn=overview_fn)
    # print(detailed_cpds)
    print("* creating detailed reports...")
    df_detailed = df[df["Compound_Id"].isin(detailed_cpds)]
    for _, rec in df_detailed.iterrows():
        cpd_id = rec["Compound_Id"]
        fn = op.join(dirname, "details", "{}.html".format(cpd_id))
        title = "{} Details".format(cpd_id)
        sim_refs = detailed_cpds[cpd_id]
        details = detailed_report(rec, sim_refs)
        write_page(details, title=title, fn=fn)

    return HTML('<a href="{}">{}</a>'.format(overview_fn, "Overview"))
