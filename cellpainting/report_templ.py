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

HTML_FILE_NAME = "report.html"
TABLE_OPTIONS = {"cellspacing": "1", "cellpadding": "1", "border": "1",
                 "align": "", "height": "60px", "summary": "Table", }  # "width": "800px",
PAGE_OPTIONS = {"icon": "icons/benzene.png"}

CSS = """
<style>
table {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border: none;
  background-color: #FFFFFF;
  text-align: left;
}
body{
  background-color: #FFFFFF;
  font-family: freesans, arial, verdana, sans-serif;
}
th {
  border-collapse: collapse;
  border-width: thin;
  border-style: solid;
  border-color: black;
  background-color: #94CAEF;
  text-align: left;
  font-weight: bold;
}
td {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  padding: 5px;
}
</style>
"""

TABLE_INTRO = """<table id="table" width="" cellspacing="1" cellpadding="1" border="1" align="center" height="60" summary="">
</tbody>"""
TABLE_EXTRO = """</tbody>
</table>"""
HTML_INTRO = """<!DOCTYPE html>
<html>
<head>
  <title>$title</title>
  <meta charset="UTF-8">
  $css
</head>
<body>
"""
t = Template(HTML_INTRO)
HTML_INTRO = t.substitute(css=CSS, title="$title")  # outwitting the formatter...
HTML_EXTRO = """
</body>
</html>"""

OVERVIEW_TABLE_HEADER = """
<tr>
    <td><b>Idx</b></td>
    <td><b>Mol</b></td>
    <td><b>Compound_Id</b></td>
    <td><b>Producer</b></td>
    <td><b>Activity<br>Flag</b></td>
    <td><b>Activity [%]</b></td>
    <td><b>Purity<br>Flag</b></td>
    <td><b>Toxic</b></td>
    <td><b>Cell<br>Fitness</b></td>
    <td><b>Number of<br>Similar References</b></td>
    <td><b>Link to<br>Detailed Report</b></td>
</tr>"""

OVERVIEW_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td>$mol_img</td>
    <td>$Compound_Id</td>
    <td>$Producer</td>
    <td>$Act_Flag</td>
    <td>$Activity</td>
    <td>$Pure_Flag</td>
    <td>$Toxic</td>
    <td>$Fitness</td>
    <td>$Num_Sim_Ref</td>
    <td>$Link</td>
</tr>"""



def page(content, title="Results", header=None, summary=None, options=PAGE_OPTIONS):
    """create a full HTML page from a list of stubs below
    options dict:
    css:     list of CSS style file paths to include.
    scripts: list of javascript library file paths to include.
    icon: path to icon image
    returns HTML page as STRING !!!"""

    # override the title if there is a title in <options>
    if "title" in options and len(options["title"]) > 2:
        title = options["title"]

    if "icon" in options and len(options["icon"]) > 2:
        icon_str = '<link rel="shortcut icon" href="{}" />'.format(options["icon"])
    else:
        icon_str = ""

    if "css" in options and options["css"]:
        css = options["css"]
        if not isinstance(css, list):
            css = [css]

        css_str = "".join(['  <link rel="stylesheet" type="text/css" href="{}">\n'.format(file_name) for file_name in css])
    else:
        # minimal inline CSS
        css_str = CSS
    if header:
        header_str = "<h2>{}</h2>".format(header)
    else:
        header_str = ""
    if summary:
        summary_str = "<p>{}</p>".format(summary)
    else:
        summary_str = ""
    if isinstance(content, list):
        content_str = "".join(content)
    else:
        content_str = content


    html_page = """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  {icon_str}
{css_str}
</head>
<body>
{header_str}
{summary_str}
{content_str}
</body>
</html>
""".format(title=title, icon_str=icon_str, css_str=css_str,
           header_str=header_str, summary_str=summary_str, content_str=content_str)

    return html_page
