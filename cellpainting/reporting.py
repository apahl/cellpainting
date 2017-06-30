#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#########
Reporting
#########

*Created on Thu Jun  8 14:40 2017 by A. Pahl*

Tools for creating HTML Reports."""

import os.path as op
import time
from string import Template
from cStringIO import StringIO as IO

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

from PIL import Image, ImageChops

from . import tools as cpt
from . import report_templ as cprt

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18


def mol_from_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles("C")
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


def write_page(page, title="Report", fn="index.html"):
    t = Template(cprt.HTML_INTRO + page + cprt.HTML_EXTRO)
    result = t.substitute(title=title)
    with open(fn, "w") as f:
        f.write(result)


def overview_report(df, cutoff=0.5):
    for _, rec in df.iterrows:
        mol = mol_from_smiles(rec["Smiles"].get("C"))



def full_report(df, dirname="report", plate=None, cutoff=0.5):
    cpt.create_dirs(op.join(dirname, "detailed"))
    print("* creating overview...")
    header = "<h2>Cell Painting Overview Report</h2>\n"
    title = ""
    if plate is not None:
        title = plate
        header += "<h3>Plate {}</h3>\n".format(plate)
    header += "<h3>({})</h3>\n".format(time.strftime("%d-%m-%Y %H:%M", time.localtime()))
    overview = header + overview_report(df, cutoff=cutoff)
    write_page(overview, title=title)



