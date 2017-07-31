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
  text-align: center;
  font-weight: bold;
  padding: 5px;
}
td {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  text-align: center;
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
  text-align: left;
  border-color:white;
  padding: 5px;
}
th.left {
    text-align: left;
}
td.left {
    text-align: left;
}
img.cpd_image {
    position: absolute;
    left: 600px;
    border: 1px solid #000000;
    padding: 10px;
}
img.logo {
    position: absolute;
    right: 0px;
    width: 300px;
}
</style>"""

JS = """
function sortOverviewTable(col_no, ascending) {
  var table, rows, switching, i, x, y, shouldSwitch;
  table = document.getElementById("overview");
  switching = true;
  /*Make a loop that will continue until
  no switching has been done:*/
  while (switching) {
    //start by saying: no switching is done:
    switching = false;
    rows = table.getElementsByTagName("TR");
    /*Loop through all table rows (except the
    first, which contains table headers):*/
    for (i = 1; i < (rows.length - 1); i++) {
      //start by saying there should be no switching:
      shouldSwitch = false;
      /*Get the two elements you want to compare,
      one from current row and one from the next:*/
      x = rows[i].getElementsByTagName("TD")[col_no];
      y = rows[i + 1].getElementsByTagName("TD")[col_no];
      //check if the two rows should switch place:
      if (ascending) {
        if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
            //if so, mark as a switch and break the loop:
            shouldSwitch = true;
            break;
        }
      } else {
        if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
            //if so, mark as a switch and break the loop:
            shouldSwitch = true;
            break;
        }
      }
    }
    if (shouldSwitch) {
      /*If a switch has been marked, make the switch
      and mark that a switch has been done:*/
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
    }
  }
}"""

lpath = op.join(op.dirname(__file__), "../res/Comas3.jpg")
logo = Image.open(lpath)
LOGO = '<img class="logo" src="data:image/jpg;base64,{}" alt="Logo"/>'
LOGO = LOGO.format(b64_img(logo))

TABLE_INTRO = """
<table id="table" width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">"""
OVERVIEW_TABLE_INTRO = """
<table id="overview" width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">"""
TABLE_EXTRO = """
</table>"""
HTML_INTRO = """
<!DOCTYPE html>
<html>
<head>
  <title>$title</title>
  <meta charset="UTF-8">
  §CSS
</head>
<body>"""
t = PreTemplate(HTML_INTRO)
HTML_INTRO = t.substitute(CSS=CSS)

OVERVIEW_HTML_INTRO = """
<!DOCTYPE html>
<html>
<head>
  <title>$title</title>
  <meta charset="UTF-8">
  §CSS
  <script>
  §JS
  </script>
</head>
<body>"""
t = PreTemplate(OVERVIEW_HTML_INTRO)
OVERVIEW_HTML_INTRO = t.substitute(CSS=CSS, JS=JS)

HTML_EXTRO = """
</body>
</html>"""

OVERVIEW_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Container_Id<br>
        <b onclick="sortOverviewTable(2, true)" title="sort ascending">&darr;</b>
        <b onclick="sortOverviewTable(2, false)" title="sort descending">&uarr;</b>
    </th>
    <th>Producer<br>
        <b onclick="sortOverviewTable(3, true)" title="sort ascending">&darr;</b>
        <b onclick="sortOverviewTable(3, false)" title="sort descending">&uarr;</b><br>
    </th>
    <th>Activity<br>Flag</th>
    <th>Induction [%]<br>
        <b onclick="sortOverviewTable(5, true)" title="sort ascending">&darr;</b>
        <b onclick="sortOverviewTable(5, false)" title="sort descending">&uarr;</b><br>
    </th>
    <th>Highest Similarity<br>to a Reference [%]<br>
        <b onclick="sortOverviewTable(6, true)" title="sort ascending">&darr;</b>
        <b onclick="sortOverviewTable(6, false)" title="sort descending">&uarr;</b><br>
    </th>
    <th>Purity<br>Flag</th>
    <th>Toxic</th>
    <th title="relative cell count in percent\nbased on the median cell count\nof the controls">
        Cell Count<br>[%Ctrl]<br>
        <b onclick="sortOverviewTable(9, true)" title="sort ascending">&darr;</b>
        <b onclick="sortOverviewTable(9, false)" title="sort descending">&uarr;</b><br>
    </th>
    <th>Link to<br>Detailed Report</th>
</tr>"""

OVERVIEW_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td><a href="http://oracle-server.mpi-dortmund.mpg.de:9944/perlbin/runjob.pl?_protocol=%7B5E6F1B34-70E1-88E4-1C51-71D239DAE9E1%7D&Compound_ID=$Compound_Id&Supplier_ID=&Alternate_ID=&IC50_view=False&Compound_level=True&Batch_level=False&__QuickRun=true">$mol_img</a></td>
    <td title="Container_Id">$Container_Id</td>
    <td title="Producer">$Producer</td>
    <td title="Activity\nFlag" bgcolor=$Col_Act_Flag>$Act_Flag</td>
    <td title="Induction [%]" bgcolor=$Col_Act>$Activity</td>
    <td title="Highest Similarity\nto a Reference [%]" bgcolor=$Col_Sim>$Max_Sim</td>
    <td title="Purity\nFlag" bgcolor=$Col_Purity>$Pure_Flag</td>
    <td title="Toxic" bgcolor=$Col_Toxic>$Toxic</td>
    <td title="Cell Count\n[%Ctrl]" bgcolor=$Col_Cell_Count>$Rel_Cell_Count</td>
    <td>$Link</td>
</tr>"""

OVERVIEW_TEMPL = """
§LOGO
"""

REF_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Container&nbsp;Id</th>
    <th>Similarity&nbsp;[%]</th>
    <th>Trivial Name</th>
    <th class="left">Known Activity</th>
</tr>"""

REF_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td><a href="../../references/details/$Container_Id.html">$mol_img</a></td>
    <td>$Container_Id</td>
    <td>$Sim_Format</td>
    <td>$Trivial_Name</td>
    <td class="left">$Known_Act</td>
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
    <td>C<br>o<br>m<br>p<br>o<br>u<br>n<br>d</td>
    <td>$Img_1_Cpd</td>
    <td>$Img_2_Cpd</td>
    <td>$Img_3_Cpd</td>
    <td>$Img_4_Cpd</td>
    <td>$Img_5_Cpd</td>
</tr>
<tr>
    <td>C<br>o<br>n<br>t<br>r<br>o<br>l</td>
    <td>$Img_1_Ctrl</td>
    <td>$Img_2_Ctrl</td>
    <td>$Img_3_Ctrl</td>
    <td>$Img_4_Ctrl</td>
    <td>$Img_5_Ctrl</td>
</tr>"""

INC_PARM_TABLE_HEADER = """
<tr>
    <th class="left">Idx</th>
    <th class="left">Increased Parameter</th>
</tr>"""

DEC_PARM_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th class="left">Decreased Parameter</th>
</tr>"""

PARM_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td class="left">$Parameter</td>
</tr>"""

DETAILS_REF_ROW = """
<tr><td class="noborder">Trivial<br>Name:</td><td class="noborder">$Trivial_Name</td></tr>
<tr><td class="noborder">Known<br>Activity:</td><td class="noborder">$Known_Act</td></tr>"""

DETAILS_TEMPL = """
§LOGO
<h1>Detailed Report</h1>
<a href="http://oracle-server.mpi-dortmund.mpg.de:9944/perlbin/runjob.pl?_protocol=%7B5E6F1B34-70E1-88E4-1C51-71D239DAE9E1%7D&Compound_ID=$Compound_Id&Supplier_ID=&Alternate_ID=&IC50_view=False&Compound_level=True&Batch_level=False&__QuickRun=true">$mol_img</a>
<h2>Compound Id $Compound_Id</h2>
§TABLE_INTRO
<tr>
    <td class="noborder" width="100px">Container Id:</td>
    <td class="noborder" width="400px">$Container_Id</td>
</tr>
<tr><td class="noborder">Producer:</td><td class="noborder">$Producer</td></tr>
<tr><td class="noborder">Induction:</td><td class="noborder">$Activity %</td></tr>
<tr><td class="noborder">Purity Flag:</td><td class="noborder">$Pure_Flag</td></tr>
<tr><td class="noborder">Cell Count:</td><td class="noborder">$Rel_Cell_Count %Ctrl</td></tr>
$Reference
§TABLE_EXTRO
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
<h3>Distribution of Active Parameters over Compartments</h3>
<p>Relative occurrence of increased and decreased parameters (compared to control)<br>
per cell compartment.</p>
$parm_hist
<br>
<br>
<h3>List of Increased Parameters (Compared to Control)</h3>
§TABLE_INTRO
§INC_PARM_TABLE_HEADER
$Inc_Parm_Table
§TABLE_EXTRO
<p></p>
<br>
<br>
<h3>List of Decreased Parameters (Compared to Control)</h3>
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
