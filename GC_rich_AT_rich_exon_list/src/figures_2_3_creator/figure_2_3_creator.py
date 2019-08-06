#!/usr/bin/env python3.5

# -*- coding: utf-8 -*-

"""
Description: The goal of this script is t ocreate the figure \
2 and 3 of the article with custom list of exons
"""


import sys
import os
import lazyparser as lp
base1 = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, base1 + "/boxplot_GC_content_and_flanking_intron_size/")
from launcher import main_1d
sys.path.insert(0, base1 + "/metaexon_figure/")
from launcher_metaexon import main as main_2a


@lp.parse(list_file="file", seddb="file", exon_type=["CCE", "ALL"],
          regulation=["up", "down"], output="dir")
def figure_creator(list_file, name_file, seddb, fasterdb, output,
                   exon_type="CCE", regulation="down", nt="S"):
    """
    Create the figure 2 and 3 with custom list of exons.

    :param list_file: (List(vtype=str)) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (List(vtype=str)) the name of each files of exons \
    given in ``list_file``
    :param seddb: (str) path to sed database
    :param exon_type: (str) the control exons
    :param regulation: (str) the resultation wanted up or down
    :param output: (str) pat were the result will be created
    :param nt: (str) the nt we want to use for the figure 1.1D
    :return:
    """
    if output[-1] != "/":
        output += "/"
    print("Creating figure 1D...")
    main_1d(list_file, name_file, seddb, exon_type, regulation, output, nt)
    print("Creation of figure 2A...")
    main_2a(list_file, name_file, "C,S,A,T,G,W,Y,R".split(","), "2A_metaexon",
            exon_type, ["#0000FF", "#00aa00"], False, output, fasterdb)


if __name__ == "__main__":
    figure_creator()