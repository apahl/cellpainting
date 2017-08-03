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

from .config import (LIMIT_ACTIVITY_H, LIMIT_ACTIVITY_L,
                     LIMIT_CELL_COUNT_H, LIMIT_CELL_COUNT_L,
                     LIMIT_SIMILARITY_H, LIMIT_SIMILARITY_L)


class PreTemplate(Template):
    delimiter = "§"


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


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

COL_WHITE = "#ffffff"
COL_GREEN = "#ccffcc"
COL_YELLOW = "#ffffcc"
COL_RED = "#ffcccc"

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
td.separate {
  border-collapse: separate;
  border-width:thin;
  border-style:solid;
  border-color:black;
  text-align: center;
  padding: 5px;
}
table.noborder {
  border-collapse: separate;
  border-width:thin;
  border-style:solid;
  border: none;
  background-color: #FFFFFF;
  text-align: left;
}
th.noborder {
  border-collapse: separate
  border-width: thin;
  border-style: solid;
  border-color: white;
  background-color: #94CAEF;
  text-align: left;
  font-weight: bold;
  padding: 5px;
}
td.noborder {
  border-collapse: separate
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
/*
  SortTable
  version 2
  7th April 2007
  Stuart Langridge, http://www.kryogenix.org/code/browser/sorttable/

  Instructions:
  Download this file
  Add class="sortable" to any table you'd like to make sortable
  Click on the headers to sort

  Thanks to many, many people for contributions and suggestions.
  Licenced as X11: http://www.kryogenix.org/code/browser/licence.html
  This basically means: do what you want with it.
*/


var stIsIE = /*@cc_on!@*/false;

sorttable = {
  init: function() {
    // quit if this function has already been called
    if (arguments.callee.done) return;
    // flag this function so we don't do the same thing twice
    arguments.callee.done = true;
    // kill the timer
    if (_timer) clearInterval(_timer);

    if (!document.createElement || !document.getElementsByTagName) return;

    sorttable.DATE_RE = /^(\\d\\d?)[\\/\\.-](\\d\\d?)[\\/\\.-]((\\d\\d)?\\d\\d)$$/;

    forEach(document.getElementsByTagName('table'), function(table) {
      if (table.className.search(/\\bsortable\\b/) != -1) {
        sorttable.makeSortable(table);
      }
    });

  },

  makeSortable: function(table) {
    if (table.getElementsByTagName('thead').length == 0) {
      // table doesn't have a tHead. Since it should have, create one and
      // put the first table row in it.
      the = document.createElement('thead');
      the.appendChild(table.rows[0]);
      table.insertBefore(the,table.firstChild);
    }
    // Safari doesn't support table.tHead, sigh
    if (table.tHead == null) table.tHead = table.getElementsByTagName('thead')[0];

    if (table.tHead.rows.length != 1) return; // can't cope with two header rows

    // Sorttable v1 put rows with a class of "sortbottom" at the bottom (as
    // "total" rows, for example). This is B&R, since what you're supposed
    // to do is put them in a tfoot. So, if there are sortbottom rows,
    // for backwards compatibility, move them to tfoot (creating it if needed).
    sortbottomrows = [];
    for (var i=0; i<table.rows.length; i++) {
      if (table.rows[i].className.search(/\\bsortbottom\\b/) != -1) {
        sortbottomrows[sortbottomrows.length] = table.rows[i];
      }
    }
    if (sortbottomrows) {
      if (table.tFoot == null) {
        // table doesn't have a tfoot. Create one.
        tfo = document.createElement('tfoot');
        table.appendChild(tfo);
      }
      for (var i=0; i<sortbottomrows.length; i++) {
        tfo.appendChild(sortbottomrows[i]);
      }
      delete sortbottomrows;
    }

    // work through each column and calculate its type
    headrow = table.tHead.rows[0].cells;
    for (var i=0; i<headrow.length; i++) {
      // manually override the type with a sorttable_type attribute
      if (!headrow[i].className.match(/\\bsorttable_nosort\\b/)) { // skip this col
        mtch = headrow[i].className.match(/\\bsorttable_([a-z0-9]+)\\b/);
        if (mtch) { override = mtch[1]; }
        if (mtch && typeof sorttable["sort_"+override] == 'function') {
          headrow[i].sorttable_sortfunction = sorttable["sort_"+override];
        } else {
          headrow[i].sorttable_sortfunction = sorttable.guessType(table,i);
        }
       // make it clickable to sort
        headrow[i].sorttable_columnindex = i;
        headrow[i].sorttable_tbody = table.tBodies[0];
        dean_addEvent(headrow[i],"click", sorttable.innerSortFunction = function(e) {

          if (this.className.search(/\\bsorttable_sorted\\b/) != -1) {
            // if we're already sorted by this column, just
            // reverse the table, which is quicker
            sorttable.reverse(this.sorttable_tbody);
            this.className = this.className.replace('sorttable_sorted',
                                                    'sorttable_sorted_reverse');
            this.removeChild(document.getElementById('sorttable_sortfwdind'));
            sortrevind = document.createElement('span');
            sortrevind.id = "sorttable_sortrevind";
            sortrevind.innerHTML = stIsIE ? '&nbsp<font face="webdings">5</font>' : '&nbsp;&#x25B4;';
            this.appendChild(sortrevind);
            return;
          }
          if (this.className.search(/\\bsorttable_sorted_reverse\\b/) != -1) {
            // if we're already sorted by this column in reverse, just
            // re-reverse the table, which is quicker
            sorttable.reverse(this.sorttable_tbody);
            this.className = this.className.replace('sorttable_sorted_reverse',
                                                    'sorttable_sorted');
            this.removeChild(document.getElementById('sorttable_sortrevind'));
            sortfwdind = document.createElement('span');
            sortfwdind.id = "sorttable_sortfwdind";
            sortfwdind.innerHTML = stIsIE ? '&nbsp<font face="webdings">6</font>' : '&nbsp;&#x25BE;';
            this.appendChild(sortfwdind);
            return;
          }

          // remove sorttable_sorted classes
          theadrow = this.parentNode;
          forEach(theadrow.childNodes, function(cell) {
            if (cell.nodeType == 1) { // an element
              cell.className = cell.className.replace('sorttable_sorted_reverse','');
              cell.className = cell.className.replace('sorttable_sorted','');
            }
          });
          sortfwdind = document.getElementById('sorttable_sortfwdind');
          if (sortfwdind) { sortfwdind.parentNode.removeChild(sortfwdind); }
          sortrevind = document.getElementById('sorttable_sortrevind');
          if (sortrevind) { sortrevind.parentNode.removeChild(sortrevind); }

          this.className += ' sorttable_sorted';
          sortfwdind = document.createElement('span');
          sortfwdind.id = "sorttable_sortfwdind";
          sortfwdind.innerHTML = stIsIE ? '&nbsp<font face="webdings">6</font>' : '&nbsp;&#x25BE;';
          this.appendChild(sortfwdind);

          // build an array to sort. This is a Schwartzian transform thing,
          // i.e., we "decorate" each row with the actual sort key,
          // sort based on the sort keys, and then put the rows back in order
          // which is a lot faster because you only do getInnerText once per row
          row_array = [];
          col = this.sorttable_columnindex;
          rows = this.sorttable_tbody.rows;
          for (var j=0; j<rows.length; j++) {
            row_array[row_array.length] = [sorttable.getInnerText(rows[j].cells[col]), rows[j]];
          }
          /* If you want a stable sort, uncomment the following line */
          //sorttable.shaker_sort(row_array, this.sorttable_sortfunction);
          /* and comment out this one */
          row_array.sort(this.sorttable_sortfunction);

          tb = this.sorttable_tbody;
          for (var j=0; j<row_array.length; j++) {
            tb.appendChild(row_array[j][1]);
          }

          delete row_array;
        });
      }
    }
  },

  guessType: function(table, column) {
    // guess the type of a column based on its first non-blank row
    sortfn = sorttable.sort_alpha;
    for (var i=0; i<table.tBodies[0].rows.length; i++) {
      text = sorttable.getInnerText(table.tBodies[0].rows[i].cells[column]);
      if (text != '') {
        if (text.match(/^-?[£$$¤]?[\\d,.]+%?$$/)) {
          return sorttable.sort_numeric;
        }
        // check for a date: dd/mm/yyyy or dd/mm/yy
        // can have / or . or - as separator
        // can be mm/dd as well
        possdate = text.match(sorttable.DATE_RE)
        if (possdate) {
          // looks like a date
          first = parseInt(possdate[1]);
          second = parseInt(possdate[2]);
          if (first > 12) {
            // definitely dd/mm
            return sorttable.sort_ddmm;
          } else if (second > 12) {
            return sorttable.sort_mmdd;
          } else {
            // looks like a date, but we can't tell which, so assume
            // that it's dd/mm (English imperialism!) and keep looking
            sortfn = sorttable.sort_ddmm;
          }
        }
      }
    }
    return sortfn;
  },

  getInnerText: function(node) {
    // gets the text we want to use for sorting for a cell.
    // strips leading and trailing whitespace.
    // this is *not* a generic getInnerText function; it's special to sorttable.
    // for example, you can override the cell text with a customkey attribute.
    // it also gets .value for <input> fields.

    if (!node) return "";

    hasInputs = (typeof node.getElementsByTagName == 'function') &&
                 node.getElementsByTagName('input').length;

    if (node.getAttribute("sorttable_customkey") != null) {
      return node.getAttribute("sorttable_customkey");
    }
    else if (typeof node.textContent != 'undefined' && !hasInputs) {
      return node.textContent.replace(/^\\s+|\\s+$$/g, '');
    }
    else if (typeof node.innerText != 'undefined' && !hasInputs) {
      return node.innerText.replace(/^\\s+|\\s+$$/g, '');
    }
    else if (typeof node.text != 'undefined' && !hasInputs) {
      return node.text.replace(/^\\s+|\\s+$$/g, '');
    }
    else {
      switch (node.nodeType) {
        case 3:
          if (node.nodeName.toLowerCase() == 'input') {
            return node.value.replace(/^\\s+|\\s+$$/g, '');
          }
        case 4:
          return node.nodeValue.replace(/^\\s+|\\s+$$/g, '');
          break;
        case 1:
        case 11:
          var innerText = '';
          for (var i = 0; i < node.childNodes.length; i++) {
            innerText += sorttable.getInnerText(node.childNodes[i]);
          }
          return innerText.replace(/^\\s+|\\s+$$/g, '');
          break;
        default:
          return '';
      }
    }
  },

  reverse: function(tbody) {
    // reverse the rows in a tbody
    newrows = [];
    for (var i=0; i<tbody.rows.length; i++) {
      newrows[newrows.length] = tbody.rows[i];
    }
    for (var i=newrows.length-1; i>=0; i--) {
       tbody.appendChild(newrows[i]);
    }
    delete newrows;
  },

  /* sort functions
     each sort function takes two parameters, a and b
     you are comparing a[0] and b[0] */
  sort_numeric: function(a,b) {
    aa = parseFloat(a[0].replace(/[^0-9.-]/g,''));
    if (isNaN(aa)) aa = 0;
    bb = parseFloat(b[0].replace(/[^0-9.-]/g,''));
    if (isNaN(bb)) bb = 0;
    return aa-bb;
  },
  sort_alpha: function(a,b) {
    if (a[0]==b[0]) return 0;
    if (a[0]<b[0]) return -1;
    return 1;
  },
  sort_ddmm: function(a,b) {
    mtch = a[0].match(sorttable.DATE_RE);
    y = mtch[3]; m = mtch[2]; d = mtch[1];
    if (m.length == 1) m = '0'+m;
    if (d.length == 1) d = '0'+d;
    dt1 = y+m+d;
    mtch = b[0].match(sorttable.DATE_RE);
    y = mtch[3]; m = mtch[2]; d = mtch[1];
    if (m.length == 1) m = '0'+m;
    if (d.length == 1) d = '0'+d;
    dt2 = y+m+d;
    if (dt1==dt2) return 0;
    if (dt1<dt2) return -1;
    return 1;
  },
  sort_mmdd: function(a,b) {
    mtch = a[0].match(sorttable.DATE_RE);
    y = mtch[3]; d = mtch[2]; m = mtch[1];
    if (m.length == 1) m = '0'+m;
    if (d.length == 1) d = '0'+d;
    dt1 = y+m+d;
    mtch = b[0].match(sorttable.DATE_RE);
    y = mtch[3]; d = mtch[2]; m = mtch[1];
    if (m.length == 1) m = '0'+m;
    if (d.length == 1) d = '0'+d;
    dt2 = y+m+d;
    if (dt1==dt2) return 0;
    if (dt1<dt2) return -1;
    return 1;
  },

  shaker_sort: function(list, comp_func) {
    // A stable sort function to allow multi-level sorting of data
    // see: http://en.wikipedia.org/wiki/Cocktail_sort
    // thanks to Joseph Nahmias
    var b = 0;
    var t = list.length - 1;
    var swap = true;

    while(swap) {
        swap = false;
        for(var i = b; i < t; ++i) {
            if ( comp_func(list[i], list[i+1]) > 0 ) {
                var q = list[i]; list[i] = list[i+1]; list[i+1] = q;
                swap = true;
            }
        } // for
        t--;

        if (!swap) break;

        for(var i = t; i > b; --i) {
            if ( comp_func(list[i], list[i-1]) < 0 ) {
                var q = list[i]; list[i] = list[i-1]; list[i-1] = q;
                swap = true;
            }
        } // for
        b++;

    } // while(swap)
  }
}

/* ******************************************************************
   Supporting functions: bundled here to avoid depending on a library
   ****************************************************************** */

// Dean Edwards/Matthias Miller/John Resig

/* for Mozilla/Opera9 */
if (document.addEventListener) {
    document.addEventListener("DOMContentLoaded", sorttable.init, false);
}

/* for Internet Explorer */
/*@cc_on @*/
/*@if (@_win32)
    document.write("<script id=__ie_onload defer src=javascript:void(0)><\\/script>");
    var script = document.getElementById("__ie_onload");
    script.onreadystatechange = function() {
        if (this.readyState == "complete") {
            sorttable.init(); // call the onload handler
        }
    };
/*@end @*/

/* for Safari */
if (/WebKit/i.test(navigator.userAgent)) { // sniff
    var _timer = setInterval(function() {
        if (/loaded|complete/.test(document.readyState)) {
            sorttable.init(); // call the onload handler
        }
    }, 10);
}

/* for other browsers */
window.onload = sorttable.init;

// written by Dean Edwards, 2005
// with input from Tino Zijdel, Matthias Miller, Diego Perini

// http://dean.edwards.name/weblog/2005/10/add-event/

function dean_addEvent(element, type, handler) {
  if (element.addEventListener) {
    element.addEventListener(type, handler, false);
  } else {
    // assign each event handler a unique ID
    if (!handler.$$$$guid) handler.$$$$guid = dean_addEvent.guid++;
    // create a hash table of event types for the element
    if (!element.events) element.events = {};
    // create a hash table of event handlers for each element/event pair
    var handlers = element.events[type];
    if (!handlers) {
      handlers = element.events[type] = {};
      // store the existing event handler (if there is one)
      if (element["on" + type]) {
        handlers[0] = element["on" + type];
      }
    }
    // store the event handler in the hash table
    handlers[handler.$$$$guid] = handler;
    // assign a global event handler to do all the work
    element["on" + type] = handleEvent;
  }
};
// a counter used to create unique IDs
dean_addEvent.guid = 1;

function removeEvent(element, type, handler) {
  if (element.removeEventListener) {
    element.removeEventListener(type, handler, false);
  } else {
    // delete the event handler from the hash table
    if (element.events && element.events[type]) {
      delete element.events[type][handler.$$$$guid];
    }
  }
};

function handleEvent(event) {
  var returnValue = true;
  // grab the event object (IE uses a global event object)
  event = event || fixEvent(((this.ownerDocument || this.document || this).parentWindow || window).event);
  // get a reference to the hash table of event handlers
  var handlers = this.events[event.type];
  // execute each event handler
  for (var i in handlers) {
    this.$$$$handleEvent = handlers[i];
    if (this.$$$$handleEvent(event) === false) {
      returnValue = false;
    }
  }
  return returnValue;
};

function fixEvent(event) {
  // add W3C standard event methods
  event.preventDefault = fixEvent.preventDefault;
  event.stopPropagation = fixEvent.stopPropagation;
  return event;
};
fixEvent.preventDefault = function() {
  this.returnValue = false;
};
fixEvent.stopPropagation = function() {
  this.cancelBubble = true;
}

// Dean's forEach: http://dean.edwards.name/base/forEach.js
/*
  forEach, version 1.0
  Copyright 2006, Dean Edwards
  License: http://www.opensource.org/licenses/mit-license.php
*/

// array-like enumeration
if (!Array.forEach) { // mozilla already supports this
  Array.forEach = function(array, block, context) {
    for (var i = 0; i < array.length; i++) {
      block.call(context, array[i], i, array);
    }
  };
}

// generic enumeration
Function.prototype.forEach = function(object, block, context) {
  for (var key in object) {
    if (typeof this.prototype[key] == "undefined") {
      block.call(context, object[key], key, object);
    }
  }
};

// character enumeration
String.forEach = function(string, block, context) {
  Array.forEach(string.split(""), function(chr, index) {
    block.call(context, chr, index, string);
  });
};

// globally resolve forEach enumeration
var forEach = function(object, block, context) {
  if (object) {
    var resolve = Object; // default
    if (object instanceof Function) {
      // functions have a "length" property
      resolve = Function;
    } else if (object.forEach instanceof Function) {
      // the object implements a custom forEach method so use that
      object.forEach(block, context);
      return;
    } else if (typeof object == "string") {
      // the object is a string
      resolve = String;
    } else if (typeof object.length == "number") {
      // the object is array-like
      resolve = Array;
    }
    resolve.forEach(object, block, context);
  }
};"""

lpath = op.join(op.dirname(__file__), "../res/Comas3.jpg")
logo = Image.open(lpath)
LOGO = '<img class="logo" src="data:image/jpg;base64,{}" alt="Logo"/>'
LOGO = LOGO.format(b64_img(logo))

TABLE_INTRO = """
<table id="table" width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">"""
OVERVIEW_TABLE_INTRO = """
<table class=sortable width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">"""
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
write(JS, "test_js.txt")
write(OVERVIEW_HTML_INTRO, "test_intro.txt")

HTML_EXTRO = """
</body>
</html>"""

OVERVIEW_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Container_Id</th>
    <th>Conc [&micro;M]</th>
    <th>Producer</th>
    <th>Activity<br>Flag</th>
    <th>Induction [%]</th>
    <th>Highest Similarity<br>to a Reference [%]</th>
    <th>Purity<br>Flag</th>
    <th>Toxic</th>
    <th title="relative cell count in percent\nbased on the median cell count\nof the controls">
        Cell Count<br>[%Ctrl]</th>
    <th>Link to<br>Detailed Report</th>
</tr>"""

OVERVIEW_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td><a href="http://oracle-server.mpi-dortmund.mpg.de:9944/perlbin/runjob.pl?_protocol=%7B5E6F1B34-70E1-88E4-1C51-71D239DAE9E1%7D&Compound_ID=$Compound_Id&Supplier_ID=&Alternate_ID=&IC50_view=False&Compound_level=True&Batch_level=False&__QuickRun=true">$mol_img</a></td>
    <td title="Container_Id">$Container_Id</td>
    <td title="Concentration">${Conc_uM}</td>
    <td title="Producer">$Producer</td>
    <td title="Activity\nFlag" bgcolor=$Col_Act_Flag>$Act_Flag</td>
    <td title="Induction [%]" bgcolor=$Col_Act>$Activity</td>
    <td title="Highest Similarity\nto a Reference [%]" bgcolor=$Col_Sim>$Max_Sim</td>
    <td title="Purity\nFlag" bgcolor=$Col_Purity>$Pure_Flag</td>
    <td title="Toxic" bgcolor=$Col_Toxic>$Toxic</td>
    <td title="Cell Count\n[%Ctrl]" bgcolor=$Col_Cell_Count>$Rel_Cell_Count</td>
    <td>$Link</td>
</tr>"""

HIGHLIGHT_LEGEND = """
<br>
<br>
<h3>Highlighting Legend</h3>
<table class=noborder id="table" width="" cellspacing="1" cellpadding="1" border="1" height="60" summary="">
<tr>
    <td class=noborder>Activity Flag:</td>
    <td bgcolor=§COL_GREEN>active</td>
    <td></td>
    <td bgcolor=§COL_RED>inactive</td>
</tr>
<tr>
    <td class=noborder>Induction [%]:</td>
    <td bgcolor=§COL_GREEN>&ge; §LIMIT_ACTIVITY_H</td>
    <td bgcolor=§COL_YELLOW>&ge; §LIMIT_ACTIVITY_L</td>
    <td bgcolor=§COL_RED>&lt; §LIMIT_ACTIVITY_L</td>
</tr>
<tr>
    <td class=noborder>Highest Similarity<br>to a Reference [%]:</td>
    <td bgcolor=§COL_GREEN>&ge; §LIMIT_SIMILARITY_H</td>
    <td bgcolor=§COL_YELLOW>&ge; §LIMIT_SIMILARITY_L</td>
    <td bgcolor=§COL_RED>&lt; §LIMIT_SIMILARITY_L</td>
</tr>
<tr>
    <td class=noborder>Purity Flag:</td>
    <td bgcolor=§COL_GREEN>Ok<br>(&ge; 80%)</td>
    <td bgcolor=§COL_YELLOW>Warn<br>(&ge; 50%)</td>
    <td bgcolor=§COL_RED>Fail<br>(&lt; 50%)</td>
</tr>
<tr>
    <td class=noborder>Toxic:</td>
    <td bgcolor=§COL_GREEN>No</td>
    <td></td>
    <td bgcolor=§COL_RED>Yes</td>
</tr>
<tr>
    <td class=noborder>Cell Count [%Ctrl]:</td>
    <td bgcolor=§COL_GREEN>&ge; §LIMIT_CELL_COUNT_H</td>
    <td bgcolor=§COL_YELLOW>&ge; §LIMIT_CELL_COUNT_L</td>
    <td bgcolor=§COL_RED>&lt; §LIMIT_CELL_COUNT_L</td>
</tr>"""
d = {
    "TABLE_INTRO": TABLE_INTRO, "TABLE_EXTRO": TABLE_EXTRO,
    "COL_GREEN": COL_GREEN, "COL_YELLOW": COL_YELLOW, "COL_RED": COL_RED,
    "LIMIT_ACTIVITY_H": LIMIT_ACTIVITY_H, "LIMIT_ACTIVITY_L": LIMIT_ACTIVITY_L,
    "LIMIT_CELL_COUNT_H": LIMIT_CELL_COUNT_H, "LIMIT_CELL_COUNT_L": LIMIT_CELL_COUNT_L,
    "LIMIT_SIMILARITY_H": LIMIT_SIMILARITY_H, "LIMIT_SIMILARITY_L": LIMIT_SIMILARITY_L
}
t = PreTemplate(HIGHLIGHT_LEGEND)
HIGHLIGHT_LEGEND = t.substitute(d)

REF_TABLE_HEADER = """
<tr>
    <th>Idx</th>
    <th>Mol</th>
    <th>Container&nbsp;Id</th>
    <th>Conc [&micro;M]</th>
    <th>Induction [%]</th>
    <th title="Profile similarity to the compound shown on top.">Similarity&nbsp;[%]</th>
    <th>Trivial Name</th>
    <th class="left">Known Activity</th>
</tr>"""

REF_TABLE_ROW = """
<tr>
    <td>$idx</td>
    <td><a href="../../references/details/$link.html">$mol_img</a></td>
    <td>$Container_Id</td>
    <td>${Conc_uM}</td>
    <td>$Activity</td>
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
<tr><td class="noborder">Conc:</td><td class="noborder">${Conc_uM} &micro;M</td></tr>
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
<br>
<br>
<h3>List of Channel Names and Compartments</h3>
§TABLE_INTRO
<tr>
    <td class="noborder" width="100px">Mito:</td>
    <td class="noborder" width="400px">Mitochondria</td>
</tr>
<tr>
    <td class="noborder">Ph_Golgi:</td>
    <td class="noborder">Golgi / Cell Membrane / Cytoskeleton</td>
</tr>
<tr>
    <td class="noborder">Syto:</td>
    <td class="noborder">Cytoplasmic RNA / Nucleoli</td>
</tr>
<tr><td class="noborder">ER:</td><td class="noborder">Endoplasmatic Reticulum</td></tr>
<tr><td class="noborder">Hoechst:</td><td class="noborder">Nuclei</td></tr>
§TABLE_EXTRO
<p>($Date)</p>"""
d = {"TABLE_INTRO": TABLE_INTRO, "TABLE_EXTRO": TABLE_EXTRO,
     "IMAGES_TABLE": IMAGES_TABLE, "LOGO": LOGO,
     "INC_PARM_TABLE_HEADER": INC_PARM_TABLE_HEADER,
     "DEC_PARM_TABLE_HEADER": DEC_PARM_TABLE_HEADER}
t = PreTemplate(DETAILS_TEMPL)
DETAILS_TEMPL = t.substitute(d)
