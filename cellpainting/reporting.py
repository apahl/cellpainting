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

import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

import numpy as np
from PIL import Image, ImageChops
import matplotlib.pyplot as plt

from . import tools as cpt
from . import report_templ as cprt
from . import processing as cpp

try:
    from . import resource_paths as cprp
except ImportError:
    from . import resource_paths_templ as cprp
    print("* Resource paths not found, stub loaded.")
    print("  Automatic loading of resources will not work,")
    print("  please have a look at resource_paths_templ.py")

from .config import (ACT_PROF_PARAMETERS, ACT_CUTOFF_PERC, ACT_CUTOFF_PERC_REF,
                     LIMIT_ACTIVITY_H, LIMIT_ACTIVITY_L,
                     LIMIT_CELL_COUNT_H, LIMIT_CELL_COUNT_L,
                     LIMIT_SIMILARITY_H, LIMIT_SIMILARITY_L)

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from IPython.core.display import HTML

try:
    from misc_tools import apl_tools
    AP_TOOLS = True
    # Library version
    VERSION = apl_tools.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))

except ImportError:
    AP_TOOLS = False
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))

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
        mol = Chem.MolFromSmiles("*")
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
    try:
        img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    except UnicodeEncodeError:
        print(Chem.MolToSmiles(mol))
        mol = Chem.MolFromSmiles("C")
        img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img.save(img_file, format='PNG')
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def b64_img(im, format="JPEG"):
    if isinstance(im, IO):
        needs_close = False
        img_file = im
    else:
        needs_close = True
        img_file = IO()
        im.save(img_file, format=format)
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    if needs_close:
        img_file.close()
    return b64


def mol_img_tag(mol, options=None):
    tag = """<img {} src="data:image/png;base64,{}" alt="Mol"/>"""
    if options is None:
        options = ""
    img_tag = tag.format(options, b64_mol(mol))
    return img_tag


def img_tag(im, format="jpeg", options=None):
    tag = """<img {} src="data:image/{};base64,{}" alt="Image"/>"""
    if options is None:
        options = ""
    b = b64_img(im, format=format)
    img_tag = tag.format(options, format.lower(), b)
    return img_tag


def load_control_images(src_dir):
    image_dir = op.join(src_dir, "images")
    ctrl_images = {}
    for ch in range(1, 6):
        im = load_image(image_dir, "H11", ch)
        ctrl_images[ch] = img_tag(im, options='style="width: 250px;"')
    return ctrl_images


def sanitize_filename(fn):
    result = fn.replace(":", "_").replace(",", "_")
    return result


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


def write_page(page, title="Report", fn="index.html", templ=cprt.HTML_INTRO):
    t = Template(templ + page + cprt.HTML_EXTRO)
    result = t.substitute(title=title)
    write(result, fn=fn)


def assign_colors(rec):
    if "Toxic" in rec:
        if rec["Toxic"]:
            rec["Col_Toxic"] = cprt.COL_RED
        else:
            rec["Col_Toxic"] = cprt.COL_GREEN
    else:
        rec["Col_Toxic"] = cprt.COL_WHITE

    if "Pure_Flag" in rec:
        if rec["Pure_Flag"] == "Ok":
            rec["Col_Purity"] = cprt.COL_GREEN
        elif rec["Pure_Flag"] == "Warn":
            rec["Col_Purity"] = cprt.COL_YELLOW
        elif rec["Pure_Flag"] == "Fail":
            rec["Col_Purity"] = cprt.COL_RED
        else:
            rec["Col_Purity"] = cprt.COL_WHITE
    else:
        rec["Col_Purity"] = cprt.COL_WHITE

    if rec["Rel_Cell_Count"] >= LIMIT_CELL_COUNT_H:
        rec["Col_Cell_Count"] = cprt.COL_GREEN
    elif rec["Rel_Cell_Count"] >= LIMIT_CELL_COUNT_L:
        rec["Col_Cell_Count"] = cprt.COL_YELLOW
    else:
        rec["Col_Cell_Count"] = cprt.COL_RED

    if rec["Activity"] >= LIMIT_ACTIVITY_H:
        rec["Col_Act"] = cprt.COL_GREEN
    elif rec["Activity"] >= LIMIT_ACTIVITY_L:
        rec["Col_Act"] = cprt.COL_YELLOW
    else:
        rec["Col_Act"] = cprt.COL_RED

    if rec["Act_Flag"] == "active":
        rec["Col_Act_Flag"] = cprt.COL_GREEN
    else:
        rec["Col_Act_Flag"] = cprt.COL_RED


def remove_colors(rec):
    for k in rec.keys():
        if k.startswith("Col_"):
            rec[k] = cprt.COL_WHITE


def overview_report(df, cutoff=LIMIT_SIMILARITY_L / 100,
                    highlight=False, mode="cpd"):
    """mode `int` displays similarities not to references but to other internal compounds
    (just displays the `Similarity` column)."""
    cpp.load_resource("SIM_REFS")
    sim_refs = cpp.SIM_REFS
    detailed_cpds = []
    if isinstance(df, cpp.DataSet):
        df = df.data
    if "ref" in mode:
        act_cutoff = ACT_CUTOFF_PERC_REF
    else:
        act_cutoff = ACT_CUTOFF_PERC
    t = Template(cprt.OVERVIEW_TABLE_HEADER)
    if "int" in mode:
        tbl_header = t.substitute(sim_entity="to another Test Compound")
    else:
        tbl_header = t.substitute(sim_entity="to a Reference")
    report = [cprt.OVERVIEW_TABLE_INTRO, tbl_header]
    row_templ = Template(cprt.OVERVIEW_TABLE_ROW)
    idx = 0
    for _, rec in df.iterrows():
        idx += 1
        well_id = rec["Well_Id"]
        mol = mol_from_smiles(rec.get("Smiles", "*"))
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx
        if "Pure_Flag" not in rec:
            rec["Pure_Flag"] = "n.d."

        rec["Act_Flag"] = "active"
        rec["Max_Sim"] = ""
        rec["Link"] = ""
        rec["Col_Sim"] = cprt.COL_WHITE
        has_details = True
        if rec["Activity"] < act_cutoff:
            has_details = False
            rec["Act_Flag"] = "inactive"
        # print(rec)
        # similar references are searched for non-toxic compounds with an activity >= LIMIT_ACTIVITY_L
        if rec["Activity"] < LIMIT_ACTIVITY_L or rec["Toxic"]:
            similars_determined = False
        else:
            similars_determined = True
        assign_colors(rec)
        convert_bool(rec, "Toxic")

        if has_details:
            detailed_cpds.append(well_id)
            details_fn = sanitize_filename(well_id)
            if rec.get("Is_Ref", False):
                plate = "references"
            else:
                plate = rec["Plate"]
            rec["Link"] = '<a href="../{}/details/{}.html">Detailed<br>Report</a>'.format(plate, details_fn)
            if similars_determined:
                if "int" in mode:
                    # similar = {"Similarity": [rec["Similarity"]]}
                    similar = pd.DataFrame({"Well_Id": [well_id], "Similarity": [rec["Similarity"]]})
                else:
                    similar = (sim_refs[sim_refs["Well_Id"] == well_id]
                               .sort_values("Similarity", ascending=False)
                               .reset_index())
                if len(similar) > 0:
                    max_sim = round(
                        similar["Similarity"][0] * 100, 1)  # first in the list has the highest similarity
                    rec["Max_Sim"] = max_sim
                    if max_sim >= LIMIT_SIMILARITY_H:
                        rec["Col_Sim"] = cprt.COL_GREEN
                    elif max_sim >= LIMIT_SIMILARITY_L:
                        rec["Col_Sim"] = cprt.COL_YELLOW
                    else:
                        rec["Col_Sim"] = cprt.COL_WHITE
                        print("ERROR: This should not happen (Max_Sim).")
                else:
                    rec["Max_Sim"] = "< {}".format(LIMIT_SIMILARITY_L)
                    rec["Col_Sim"] = cprt.COL_RED

        if not highlight:
            # remove all coloring again:
            remove_colors(rec)
        report.append(row_templ.substitute(rec))
    report.append(cprt.TABLE_EXTRO)
    return "\n".join(report), detailed_cpds


def sim_ref_table(similar):
    cpp.load_resource("REFERENCES")
    df_refs = cpp.REFERENCES
    table = [cprt.TABLE_INTRO, cprt.REF_TABLE_HEADER]
    templ = Template(cprt.REF_TABLE_ROW)
    for idx, rec in similar.iterrows():
        rec = rec.to_dict()
        ref_id = rec["Ref_Id"]
        ref_data = df_refs[df_refs["Well_Id"] == ref_id]
        if len(ref_data) == 0:
            raise ValueError("BUG: ref_data should not be empty.")
        ref_data = ref_data.copy()
        ref_data = ref_data.fillna("&mdash;")
        rec.update(ref_data.to_dict("records")[0])
        mol = mol_from_smiles(rec.get("Smiles", "*"))
        rec["Sim_Format"] = "{:.1f}".format(rec["Similarity"] * 100)
        rec["Tan_Format"] = "{:.1f}".format(rec["Tanimoto"] * 100)
        if rec["Tan_Format"] == np.nan:
            rec["Tan_Format"] = "&mdash;"
        rec["mol_img"] = mol_img_tag(mol)
        rec["idx"] = idx + 1

        link = sanitize_filename(rec["Well_Id"])
        rec["link"] = link
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
            "Parameter": p[7:]
        }
        row = templ.substitute(rec)
        table.append(row)
    return "\n".join(table), changed


def parm_stats(parameters):
    result = []
    channels = ["_Mito", "_Ph_golgi", "_Syto", "_ER", "Hoechst"]
    for ch in channels:
        cnt = len([p for p in parameters if ch in p])
        result.append(cnt)
    return result


def parm_hist(increased, decreased):
    labels = [
        "Mito",
        "Golgi / Membrane",
        "RNA / Nucleoli",
        "ER",
        "Nuclei"
    ]

    inc_max = max(increased)
    dec_max = max(decreased)
    max_total = max([inc_max, dec_max])
    if max_total == 0:
        result = "No compartment-specific parameters were changed."
        return result
    inc_norm = [v / max_total for v in increased]
    dec_norm = [v / max_total for v in decreased]

    n_groups = 5
    dpi = 96
    # plt.rcParams['axes.titlesize'] = 25
    plt.style.use("seaborn-white")
    plt.style.use("seaborn-pastel")
    plt.style.use("seaborn-talk")
    plt.rcParams['axes.labelsize'] = 25
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    plt.rcParams['legend.fontsize'] = 20
    size = (1500, 1000)
    figsize = (size[0] / dpi, size[1] / dpi)
    fig, ax = plt.subplots(figsize=figsize)
    index = np.arange(n_groups)
    bar_width = 0.25
    plt.bar(index, inc_norm, bar_width,
            color='#94caef',
            label='Inc')
    plt.bar(index + bar_width, dec_norm, bar_width,
            color='#ffdd1a',
            label='Dec')

    plt.xlabel('Cell Compartment')
    plt.ylabel('rel. Occurrence')
    plt.xticks(index + bar_width / 2, labels, rotation=45)
    plt.legend()
    plt.tight_layout()
    img_file = IO()
    plt.savefig(img_file, bbox_inches='tight', format="jpg")
    result = img_tag(img_file, format="jpg", options='style="width: 800px;"')
    img_file.close()
    # important, otherwise the plots will accumulate and fill up memory:
    plt.close()
    return result


def show_images(plate_quad, well):
    """For interactive viewing in the notebook."""
    date = cprp.DATES[plate_quad]
    src_dir = cprp.src_path.format(date, plate_quad)
    ctrl_images = load_control_images(src_dir)
    image_dir = op.join(src_dir, "images")
    templ_dict = {}
    for ch in range(1, 6):
        im = load_image(image_dir, well, ch)
        templ_dict["Img_{}_Cpd".format(ch)] = img_tag(im, options='style="width: 250px;"')
        templ_dict["Img_{}_Ctrl".format(ch)] = ctrl_images[ch]
    tbody_templ = Template(cprt.IMAGES_TABLE)
    table = cprt.TABLE_INTRO + tbody_templ.substitute(templ_dict) + cprt.HTML_EXTRO
    return HTML(table)


def detailed_report(rec, src_dir, ctrl_images):
    cpp.load_resource("SIM_REFS")
    sim_refs = cpp.SIM_REFS
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
    image_dir = op.join(src_dir, "images")
    well_id = rec["Well_Id"]
    act_prof = rec["Act_Profile"]
    inc_parm, changed = changed_parameters_table(act_prof, "2")
    increased = parm_stats(changed)
    dec_parm, changed = changed_parameters_table(act_prof, "0")
    decreased = parm_stats(changed)
    mol = mol_from_smiles(rec.get("Smiles", "*"))
    if "Pure_Flag" not in rec:
        rec["Pure_Flag"] = "n.d."

    templ_dict = rec.copy()
    templ_dict["Date"] = date
    templ_dict["mol_img"] = mol_img_tag(mol, options='class="cpd_image"')
    templ_dict["Inc_Parm_Table"] = inc_parm
    templ_dict["Dec_Parm_Table"] = dec_parm
    templ_dict["parm_hist"] = parm_hist(increased, decreased)
    if "Known_Act" in templ_dict:
        if templ_dict["Trivial_Name"] == np.nan or templ_dict["Trivial_Name"] == "":
            templ_dict["Trivial_Name"] = "&mdash;"
        if templ_dict["Trivial_Name"] == np.nan or templ_dict["Known_Act"] == "":
            templ_dict["Known_Act"] = "&mdash;"
        t = Template(cprt.DETAILS_REF_ROW)
        templ_dict["Reference"] = t.substitute(templ_dict)
    else:
        templ_dict["Reference"] = ""
    well = rec["Metadata_Well"]
    for ch in range(1, 6):
        im = load_image(image_dir, well, ch)
        templ_dict["Img_{}_Cpd".format(ch)] = img_tag(im, options='style="width: 250px;"')
        templ_dict["Img_{}_Ctrl".format(ch)] = ctrl_images[ch]
    if rec["Activity"] < LIMIT_ACTIVITY_L:
        templ_dict["Ref_Table"] = "Because of low induction (&lt; {}%), no similarity was determined.".format(LIMIT_ACTIVITY_L)
    elif rec["Rel_Cell_Count"] < LIMIT_CELL_COUNT_L:
        templ_dict["Ref_Table"] = "Because of compound toxicity, no similarity was determined."
    else:
        similar = (sim_refs[sim_refs["Well_Id"] == well_id]
                   .sort_values("Similarity", ascending=False)
                   .reset_index())
        if len(similar) > 0:
            ref_tbl = sim_ref_table(similar)
            templ_dict["Ref_Table"] = ref_tbl
        else:
            templ_dict["Ref_Table"] = "No similar references found."

    t = Template(cprt.DETAILS_TEMPL)
    report = t.substitute(templ_dict)
    return report


def full_report(df, src_dir, report_name="report", plate=None,
                cutoff=0.6, act_cutoff=ACT_CUTOFF_PERC, highlight=False, mode="cpd"):
    overview_fn = op.join(report_name, "index.html")
    date = time.strftime("%d-%m-%Y %H:%M", time.localtime())
    cpt.create_dirs(op.join(report_name, "details"))
    if isinstance(df, cpp.DataSet):
        df = df.data
    print("* creating overview...")
    header = "{}\n<h2>Cell Painting Overview Report</h2>\n".format(cprt.LOGO)
    title = "Overview"
    if plate is not None:
        title = plate
        header += "<h3>Plate {}</h3>\n".format(plate)
    header += "<p>({})</p>\n".format(date)
    if highlight:
        highlight_legend = cprt.HIGHLIGHT_LEGEND
    else:
        highlight_legend = ""
    overview, detailed_cpds = overview_report(
        df, cutoff=cutoff, highlight=highlight, mode=mode)
    overview = header + overview + highlight_legend
    write_page(overview, title=title, fn=overview_fn, templ=cprt.OVERVIEW_HTML_INTRO)
    # print(detailed_cpds)
    print("* creating detailed reports...")
    print("  * loading control images...")
    ctrl_images = load_control_images(src_dir)
    print("  * writing individual reports...")
    df_detailed = df[df["Well_Id"].isin(detailed_cpds)]
    for _, rec in df_detailed.iterrows():
        well_id = rec["Well_Id"]
        fn = op.join(report_name, "details", "{}.html".format(sanitize_filename(well_id)))
        title = "{} Details".format(well_id)
        # similar = detailed_cpds[well_id]
        details = detailed_report(rec, src_dir, ctrl_images)
        write_page(details, title=title, fn=fn, templ=cprt.DETAILS_HTML_INTRO)

    print("* done.")
    return HTML('<a href="{}">{}</a>'.format(overview_fn, "Overview"))
