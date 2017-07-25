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

COL_WHITE = "#ffffff"
COL_GREEN = "#77ff33"
COL_YELLOW = "#ffff99"
COL_RED = "#ff6666"


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


def load_image(path, well, channel):
    image_fn = "{}/{}_w{}.jpg".format(path, well, channel)
    im = Image.open(image_fn)
    return im


def b64_mol(mol, size=300):
    img_file = IO()
    img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img.save(img_file, format='PNG')
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def b64_img(im, format="JPEG"):
    img_file = IO()
    im.save(img_file, format=format)
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def mol_img_tag(mol, style=None):
    if style is None:
        style = ""
    else:
        style = """style='{}' """.format(style)
    img_tag = '<img {}src="data:image/png;base64,{}" alt="Mol"/>'.format(style, b64_mol(mol))
    return img_tag


def img_tag(im, format="jpeg", style=None):
    if style is None:
        style = ""
    else:
        style = """style='{}' """.format(style)
    b = b64_img(im, format=format)
    img_tag = '<img {}src="data:image/{};base64,{}" alt="Cell"/>'.format(style, format.lower(), b)
    return img_tag


def load_control_images(src_dir):
    image_dir = op.join(src_dir, "images")
    ctrl_images = {}
    for ch in range(1, 6):
        im = load_image(image_dir, "H11", ch)
        ctrl_images[ch] = img_tag(im, style='width: 250px;')
    return ctrl_images


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


def write_page(page, title="Report", fn="index.html"):
    t = Template(cprt.HTML_INTRO + page + cprt.HTML_EXTRO)
    result = t.substitute(title=title)
    write(result, fn=fn)


def assign_colors(rec):
    if "Toxic" in rec:
        if rec["Toxic"]:
            rec["Col_Toxic"] = COL_RED
        else:
            rec["Col_Toxic"] = COL_GREEN
    else:
        rec["Col_Toxic"] = COL_WHITE

    if "Pure_Flag" in rec:
        if rec["Pure_Flag"] == "Ok":
            rec["Col_Purity"] = COL_GREEN
        elif rec["Pure_Flag"] == "Warn":
            rec["Col_Purity"] = COL_YELLOW
        elif rec["Pure_Flag"] == "Fail":
            rec["Col_Purity"] = COL_RED
        else:
            rec["Col_Purity"] = COL_WHITE
    else:
        rec["Col_Purity"] = COL_WHITE

    if rec["Fitness"] >= 75:
        rec["Col_Fitness"] = COL_GREEN
    elif rec["Fitness"] >= 55:
        rec["Col_Fitness"] = COL_YELLOW
    else:
        rec["Col_Fitness"] = COL_RED

    if rec["Activity"] < ACT_CUTOFF_PERC:
        rec["Col_Active"] = COL_RED
    else:
        rec["Col_Active"] = COL_GREEN


def remove_colors(rec):
    for k in rec.keys():
        if k.startswith("Col_"):
            rec[k] = COL_WHITE


def overview_report(df, df_refs=None, cutoff=0.6, detailed_cpds=None, highlight=False):
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
        if "Pure_Flag" not in rec:
            rec["Pure_Flag"] = "n.d."
        assign_colors(rec)
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
                rec["Col_Sim"] = COL_RED
            else:
                # rec["Num_Sim_Ref"] = str(sim_refs.shape[0])
                max_sim = sim_refs["Similarity"].max() * 100
                rec["Max_Sim"] = str(max_sim)
                if max_sim >= 80:
                    rec["Col_Sim"] = COL_GREEN
                elif max_sim >= cutoff * 100:
                    rec["Col_Sim"] = COL_YELLOW
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
        if not highlight:
            # remove all coloring again:
            remove_colors(rec)
        report.append(row_templ.substitute(rec))
    report.append(cprt.TABLE_EXTRO)
    return "\n".join(report)


def sim_ref_table(sim_refs):
    table = [cprt.TABLE_INTRO, cprt.REF_TABLE_HEADER]
    templ = Template(cprt.REF_TABLE_ROW)
    for idx, sim_ref in enumerate(sim_refs, 1):
        rec = sim_ref.copy()
        mol = mol_from_smiles(rec.get("Smiles", "C"))
        rec["Similarity"] *= 100
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx
        row = templ.substitute(rec)
        table.append(row)
    table.append(cprt.TABLE_EXTRO)
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


def detailed_report(rec, sim_refs, src_dir, ctrl_images):
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
    image_dir = op.join(src_dir, "images")
    act_prof = rec["Act_Profile"]
    inc_parm = changed_parameters_table(act_prof, "2")
    dec_parm = changed_parameters_table(act_prof, "0")
    mol = mol_from_smiles(rec.get("Smiles", "C"))
    if "Pure_Flag" not in rec:
        rec["Pure_Flag"] = "n.d."

    templ_dict = {
        "Compound_Id": rec["Compound_Id"],
        "Date": date,
        "Producer": rec["Producer"],
        "Activity": rec["Activity"],
        "Fitness": rec["Fitness"],
        "Pure_Flag": rec["Pure_Flag"],
        "mol_img": mol_img_tag(mol, style='border:1px solid #000000; padding: 10px;'),
        "Inc_Parm_Table": inc_parm,
        "Dec_Parm_Table": dec_parm,
    }
    well = rec["Metadata_Well"]
    for ch in range(1, 6):
        im = load_image(image_dir, well, ch)
        templ_dict["Img_{}_Cpd".format(ch)] = img_tag(im, style='width: 250px;')
        templ_dict["Img_{}_Ctrl".format(ch)] = ctrl_images[ch]
    if len(sim_refs) > 0:
        ref_tbl = sim_ref_table(sim_refs)
        templ_dict["Ref_Table"] = ref_tbl
    else:
        templ_dict["Ref_Table"] = "No similar References found."
    t = Template(cprt.DETAILS_TEMPL)

    report = t.substitute(templ_dict)
    return report


def full_report(df, src_dir, df_refs=None, report_name="report", plate=None, cutoff=0.6, highlight=False):
    overview_fn = op.join(report_name, "index.html")
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
    cpt.create_dirs(op.join(report_name, "details"))
    if isinstance(df, cpp.DataSet):
        df = df.data
    df_refs = cpt.check_df(df_refs, REFERENCES)
    detailed_cpds = {}
    print("* creating overview...")
    header = "{}\n<h2>Cell Painting Overview Report</h2>\n".format(cprt.LOGO)
    title = "Overview"
    if plate is not None:
        title = plate
        header += "<h3>Plate {}</h3>\n".format(plate)
    header += "<p>({})</p>\n".format(date)
    overview = header + overview_report(df, df_refs, cutoff=cutoff,
                                        detailed_cpds=detailed_cpds, highlight=highlight)
    write_page(overview, title=title, fn=overview_fn)
    # print(detailed_cpds)
    print("* creating detailed reports...")
    print("  * loading Control images...")
    ctrl_images = load_control_images(src_dir)
    print("  * writing individual reports...")
    df_detailed = df[df["Compound_Id"].isin(detailed_cpds)]
    for _, rec in df_detailed.iterrows():
        cpd_id = rec["Compound_Id"]
        fn = op.join(report_name, "details", "{}.html".format(cpd_id))
        title = "{} Details".format(cpd_id)
        sim_refs = detailed_cpds[cpd_id]
        details = detailed_report(rec, sim_refs, src_dir, ctrl_images)
        write_page(details, title=title, fn=fn)

    return HTML('<a href="{}">{}</a>'.format(overview_fn, "Overview"))
