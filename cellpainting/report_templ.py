#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
################
Report Templates
################

*Created on Thu Jun  8 14:40 2017 by A. Pahl*

Templates for the HTML Reports.

from string import Template
In [2]: report_all = Template("$name, this is a $what for Templates.")
In [3]: tmap = {"name": "Axel", "what": "Test"}
In [4]: report_all.substitute(tmap)
Out[4]: 'Axel, this is a Test for Templates.'
"""

from string import Template


class PreTemplate(Template):
    delimiter = "§"


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
</style>"""

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
    <td title="Purity\nFlag">$Pure_Flag</td>
    <td title="Toxic">$Toxic</td>
    <td title="Cell\nViability [%]">$Fitness</td>
    <td title="Highest Similarity\nto a Reference [%]">$Max_Sim</td>
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
<p style="float: right;">$mol_img</p>
<h1>Detailed Report</h1>
<h2>Compound Id $Compound_Id</h2>
<p>Producer: $Producer</p>
<p>Activity: $Activity %</p>
<br>
<br>
<h3>Similar References</h3>
§TABLE_INTRO
§REF_TABLE_HEADER
$Ref_Table
§TABLE_EXTRO
<p></p>
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
     "REF_TABLE_HEADER": REF_TABLE_HEADER, "INC_PARM_TABLE_HEADER": INC_PARM_TABLE_HEADER,
     "DEC_PARM_TABLE_HEADER": DEC_PARM_TABLE_HEADER}
t = PreTemplate(DETAILS_TEMPL)
DETAILS_TEMPL = t.substitute(d)

DETAILS_TEMPL_NO_SIM_REFS = """
<p style="float: right;">$mol_img</p>
<h1>Detailed Report</h1>
<h2>Compound Id $Compound_Id</h2>
<p>Producer: $Producer</p>
<p>Activity: $Activity %</p>
<br>
<br>
<h3>Similar References</h3>
No Similar References were found.
<p></p>
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
     "INC_PARM_TABLE_HEADER": INC_PARM_TABLE_HEADER,
     "DEC_PARM_TABLE_HEADER": DEC_PARM_TABLE_HEADER}
t = PreTemplate(DETAILS_TEMPL_NO_SIM_REFS)
DETAILS_TEMPL_NO_SIM_REFS = t.substitute(d)
