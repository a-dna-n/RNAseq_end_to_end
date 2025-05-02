#!/usr/bin/env python

"""
Functions to retrieve RNA-seq data and metadata, launch alignments, counts, reassembly, etc. 

GEOparse and pySRAdb have excellent features, but neither knows about NCBI human recounts, and neither lets us choose which supplementary files to download depending on their contents. GEOparse only uses soft-formatted files, which are problematic in studies with human and mouse samples.

GSE154783: SRP272683
GSE228268: SRP429497  = c. sabaeus
"""
import sys
import os
from modules.tools import *


test_executables("free gunzip hostname pysradb samtools STAR wget".split(), exit_on_error = False)
test_libraries("numpy pandas pysradb metapub cyclopts".split(), exit_on_error = True) # GEOparse 

import cyclopts

from modules import metadata, workflow, plots, pipeline

cyc_app = cyclopts.App(help = "Functions for cyno RNA-seq analysis.")
cyc_plots = cyclopts.Group.create_ordered("plots")

cyc_app.update(metadata.cyc_app)
cyc_app.update(workflow.cyc_app)
cyc_app.update(pipeline.cyc_app)
cyc_app.update(plots.cyc_app)

"""


import re
import glob
from textwrap import dedent
from collections import OrderedDict, Counter
from types import SimpleNamespace
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
#import GEOparse
#from modules.arg_def import *
from typing import Literal, Optional # TypeAlias, Union
from modules.constants import *
#from lxml import html
#import requests
constants = define_constants()

from cyclopts import App as cyc_App, Group as cyc_Group, Parameter
from dataclasses import dataclass, KW_ONLY, field

cyc_app = cyc_App(help = "Functions for RNA-seq data and metadata.")
# function groups
cyc_groups = {}
cyc_groups["data"] = cyc_Group.create_ordered("data")
buildrefs = cyc_Group.create_ordered("buildrefs")
align = cyc_Group.create_ordered("align")
counts = cyc_Group.create_ordered("counts")
assembly = cyc_Group.create_ordered("assembly")
bamfiles = cyc_Group.create_ordered("bamfiles")
cyc_groups["utils"] = cyc_Group.create_ordered("utils")



"""


if __name__ == "__main__":
    #constants.star_indexes = find_star_indexes()
    cyc_app()
