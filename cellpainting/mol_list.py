#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
########
Mol_List
########

*Created on Thu Jun  22 11:00 2017 by A. Pahl*

Extending the rdkit_ipynb_tools.Mol_List for CellPainting.
"""

from rdkit_ipynb_tools import tools


class Mol_List(tools.Mol_List):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def _pass_properties(self, new_list):
            new_list.order = self.order
            new_list.ia = self.ia
            new_list.plot_tool = self.plot_tool


    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        try:
            new_list = type(self)(result)

            # pass on properties
            self._pass_properties(new_list)
            return new_list
        except TypeError:
            return result


    def new(self, *args):
        new_list = type(self)(*args)
        # pass on properties
        self._pass_properties(new_list)
        return new_list


    def num_of_similars_in_refs(self, ref_file="XX"):
        """Find the number of similar references for the quick report.
        Adds the number to the molecules."""
        pass
