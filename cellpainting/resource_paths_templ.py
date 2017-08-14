#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############
Resource Paths
##############

*Created on Sun Aug 6, 2017 18:30 by A. Pahl*

This is a template file. Adapt to your needs and rename to resource_paths.py
"""

# tab-delim., gzipped
smiles_path = "xxx/smiles_b64.tsv.gz"
smiles_cols = ['Compound_Id', "Smiles"]

# container based data, like purity; tab-delim., gzipped
data_path = "xxx/data_b64.tsv.gz"
data_cols = ['Container_Id', "Pure_Flag", "Avail"]

# annotations of known references
# should contain Compound_Id, Trivial_Name, Known_Act; tab-delim.
annotations_path = "xxx/known_act.tsv"

# the prepared file of references
references_path = "xxx/references_act_prof.tsv"

# storage of Container_Id and their highest Sim reference
sim_refs_path = "xxx/sim_refs.pkl"

sim_refs_ext_path = "/home/pahl/comas/notebooks/projects/painting/sim_refs_ext.pkl"

# collection of all data
datastore_path = "/home/pahl/comas/notebooks/projects/painting/cp_datastore.tsv"
layouts_path = "/home/pahl/comas/projects/painting/layouts/layouts.tsv"


# just used for viewing images in the Jupyter notebook
# date and plate_quad are used for formatting
src_path = "xxx/{}-{}"


DATES = {
}

CONFS = {
}
