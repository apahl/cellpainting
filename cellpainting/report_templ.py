#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
################
Report Templates
################

*Created on Thu Jun  8 14:40 2017 by A. Pahl*

Templates for the HTML Reports."""

import base64
import os.path as op
from string import Template
from io import BytesIO as IO

from PIL import Image


class PreTemplate(Template):
    delimiter = "§"


def b64_img(im, format="JPEG"):
    img_file = IO()
    im.save(img_file, format=format)
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


HTML_FILE_NAME = "report.html"
TABLE_OPTIONS = {"cellspacing": "1", "cellpadding": "1", "border": "1",
                 "align": "", "height": "60px", "summary": "Table", }  # "width": "800px",
PAGE_OPTIONS = {"icon": "icons/benzene.png"}

CSS = """
<style>
body {
  text-align: left;
  background-color: #FFFFFF;
  font-family: freesans, arial, verdana, sans-serif;
}
table {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border: none;
  background-color: #FFFFFF;
  text-align: left;
}
th {
  border-collapse: collapse;
  border-width: thin;
  border-style: solid;
  border-color: black;
  background-color: #94CAEF;
  text-align: left;
  font-weight: bold;
  padding: 5px;
}
td {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  padding: 5px;
}
th.noborder {
  border-collapse: collapse;
  border-width: thin;
  border-style: solid;
  border-color: white;
  background-color: #94CAEF;
  text-align: left;
  font-weight: bold;
  padding: 5px;
}
td.noborder {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:white;
  padding: 5px;
}
</style>"""

lpath = op.join(op.dirname(__file__), "../res/Comas3.png")
logo = Image.open(lpath)
LOGO = '<img style="float: right; width: 300px;" src="data:image/png;base64,{}" alt="Cell"/>'
LOGO = LOGO.format(b64_img(logo, format="PNG"))

TABLE_INTRO = """
<table id="table" width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">"""
TABLE_EXTRO = """
</table>"""
HTML_INTRO = """
<!DOCTYPE html>
<html>
<head>
  <title>$title</title>
  <meta charset="UTF-8">
  §css
</head>
<body>"""
t = PreTemplate(HTML_INTRO)
HTML_INTRO = t.substitute(css=CSS)

HTML_EXTRO = """
</body>
</html>"""

OVERVIEW_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Compound_Id</th>
    <th>Producer</th>
    <th>Activity<br>Flag</th>
    <th>Activity [%]</th>
    <th>Purity<br>Flag</th>
    <th>Toxic</th>
    <th>Cell<br>Viability [%]</th>
    <th>Highest Similarity<br>to a Reference [%]</th>
    <th>Link to<br>Detailed Report</th>
</tr>"""

OVERVIEW_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td>$mol_img</td>
    <td title="Compound_Id">$Compound_Id</td>
    <td title="Producer">$Producer</td>
    <td title="Activity\nFlag">$Act_Flag</td>
    <td title="Activity [%]">$Activity</td>
    <td title="Purity\nFlag" bgcolor=$Col_Purity>$Pure_Flag</td>
    <td title="Toxic" bgcolor=$Col_Toxic>$Toxic</td>
    <td title="Cell\nViability [%]" bgcolor=$Col_Fitness>$Fitness</td>
    <td title="Highest Similarity\nto a Reference [%]" bgcolor=$Col_Sim>$Max_Sim</td>
    <td>$Link</td>
</tr>"""

REF_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Compound_Id</th>
    <th>Similarity [%]</th>
    <th>Trivial Name</th>
    <th>Known Activity</th>
</tr>"""

REF_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td>$mol_img</td>
    <td>$Compound_Id</td>
    <td>$Similarity</td>
    <td>$Trivial_Name</td>
    <td>$Known_Act</td>
</tr>"""

IMAGES_TABLE = """
<tr>
    <th></th>
    <th>Mitochondria</th>
    <th>Golgi / Cell <br>Membrane / Cytoskeleton</th>
    <th>Cytoplasmic RNA / <br>Nucleoli</th>
    <th>Endoplasmatic<br>Reticulum</th>
    <th>Nuclei</th>
</tr>
<tr>
    <td align="center">C<br>o<br>m<br>p<br>o<br>u<br>n<br>d</td>
    <td>$Img_1_Cpd</td>
    <td>$Img_2_Cpd</td>
    <td>$Img_3_Cpd</td>
    <td>$Img_4_Cpd</td>
    <td>$Img_5_Cpd</td>
</tr>
<tr>
    <td align="center">C<br>o<br>n<br>t<br>r<br>o<br>l</td>
    <td>$Img_1_Ctrl</td>
    <td>$Img_2_Ctrl</td>
    <td>$Img_3_Ctrl</td>
    <td>$Img_4_Ctrl</td>
    <td>$Img_5_Ctrl</td>
</tr>"""

INC_PARM_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Increased Parameter</th>
</tr>"""

DEC_PARM_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Decreased Parameter</th>
</tr>"""

PARM_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td>$Parameter</td>
</tr>"""

DETAILS_TEMPL = """
§LOGO
<div>
<div style="float: left; width: 400px;">
<h1>Detailed Report</h1>
<h2>Compound Id $Compound_Id</h2>
§TABLE_INTRO
<tr><td class="noborder">Producer:</td><td class="noborder">$Producer</td></tr>
<tr><td class="noborder">Activity:</td><td class="noborder">$Activity %</td></tr>
<tr><td class="noborder">Purity Flag:</td><td class="noborder">$Pure_Flag</td></tr>
<tr><td class="noborder">Cell Viability:</td><td class="noborder">$Fitness %</td></tr>
§TABLE_EXTRO
</div>
<p><br><br>$mol_img</p>
</div>
<br>
<br>
<h3>Similar References</h3>
$Ref_Table
<p></p>
<br>
<br>
<h3>Images</h3>
Sample images from site 5.
§TABLE_INTRO
§IMAGES_TABLE
§TABLE_EXTRO
<br>
<br>
<h3>Parameters Increased Compared to Control</h3>
§TABLE_INTRO
§INC_PARM_TABLE_HEADER
$Inc_Parm_Table
§TABLE_EXTRO
<p></p>
<br>
<br>
<h3>Parameters Decreased Compared to Control</h3>
§TABLE_INTRO
§DEC_PARM_TABLE_HEADER
$Dec_Parm_Table
§TABLE_EXTRO
<p>($Date)</p>"""
d = {"TABLE_INTRO": TABLE_INTRO, "TABLE_EXTRO": TABLE_EXTRO,
     "IMAGES_TABLE": IMAGES_TABLE, "LOGO": LOGO,
     "INC_PARM_TABLE_HEADER": INC_PARM_TABLE_HEADER,
     "DEC_PARM_TABLE_HEADER": DEC_PARM_TABLE_HEADER}
t = PreTemplate(DETAILS_TEMPL)
DETAILS_TEMPL = t.substitute(d)
