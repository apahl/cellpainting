#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####################
DataFrame Web Browser
####################

*Created on Wed Aug 30 9:00 2017 by A. Pahl*

View Pandas DataFrames in the Browser.
Code is based on bluenote10's NimData html browser
(https://github.com/bluenote10/NimData/blob/master/src/nimdata/html.nim)."""

import os
from string import Template
import webbrowser

LIB_LOCATION = "local"
SHOW_WARN = True

PANDAS_TABLE_NET = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>$title</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.13/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://bootswatch.com/cosmo/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">

    <style>
      .custom {
        font-size: 12px;
      }
      .bottomcustom, .topcustom {
        font-size: 10px;
      }
    </style>
  </head>
  <body>

    <div class="container">
      <div class="page-header">
          <h3>$title</h3>
      </div>
      <table id="maintable" class="display custom" cellspacing="0" width="100%">
        $table
      </table>
    </div>

    <!-- jQuery -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
    <!-- Bootstrap -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
    <!-- Data Table -->
    <script src="https://cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js"></script>
    <script>
        $$(document).ready(function() {
            $$('#maintable').DataTable({
                "pageLength": 25,
                "dom": '<"topcustom"lfr>t<"bottomcustom"ip>'
            });
        });
    </script>
  </body>
</html>
"""

PANDAS_TABLE_LOCAL = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>$title</title>
    <link rel="stylesheet" href="lib/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap1.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap2.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap-theme.min.css">

    <style>
      .custom {
        font-size: 12px;
      }
      .bottomcustom, .topcustom {
        font-size: 10px;
      }
    </style>
  </head>
  <body>

    <div class="container">
      <div class="page-header">
          <h3>$title</h3>
      </div>
      <table id="maintable" class="display custom" cellspacing="0" width="100%">
        $table
      </table>
    </div>

    <!-- jQuery -->
    <script src="lib/jquery.min.js"></script>
    <!-- Bootstrap -->
    <script src="lib/bootstrap.min.js"></script>
    <!-- Data Table -->
    <script src="lib/jquery.dataTables.min.js"></script>
    <script>
        $$(document).ready(function() {
            $$('#maintable').DataTable({
                "pageLength": 25,
                "dom": '<"topcustom"lfr>t<"bottomcustom"ip>'
            });
        });
    </script>
  </body>
</html>
"""


def browse(df, title="DataFrame", drop=None, keep=None, fn="tmp.html"):
    global SHOW_WARN
    if "local" in LIB_LOCATION.lower() and os.access("lib/bootstrap.min.js", os.R_OK):
        pandas_tbl = PANDAS_TABLE_LOCAL
    else:
        pandas_tbl = PANDAS_TABLE_NET
        if SHOW_WARN:
            SHOW_WARN = False
            print("* using online libs for dataframe browsing...")
    if drop is not None:
        df = df.drop(drop, axis=1)
    if keep is not None:
        df = df[keep]
    tbl = df.to_html()
    tbl_list = tbl.split("\n")
    tbl_list = tbl_list[1:-1]
    tbl = "\n".join(tbl_list)
    t = Template(pandas_tbl)
    html = t.substitute(title=title, table=tbl)
    with open(fn, "w") as f:
        f.write(html)
    webbrowser.open_new_tab(fn)
