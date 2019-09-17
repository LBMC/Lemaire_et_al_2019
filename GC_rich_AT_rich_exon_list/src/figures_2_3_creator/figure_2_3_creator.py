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
base2 = os.path.dirname(os.path.dirname(base1))
sys.path.insert(0, base1 + "/boxplot_GC_content_and_flanking_intron_size/")
from launcher import main_1d
sys.path.insert(0, base1 + "/metaexon_figure/")
from launcher_metaexon import main as main_2a
sys.path.insert(0, base2 + "/Figure_ESA/src/")
from heatmap_maker import main_2bc
sys.path.insert(0, base1 + "/minimum_free_energy/")
from mfe_figure_creator import main_2d
sys.path.insert(0, base1 + "/make_control_files_bp_ppt/")
from bp_ppt_figure_creator import main_2efg
sys.path.insert(0, base1 + "/stretch_calculator/")
from figure_creator_stretch import main_2g
sys.path.insert(0, base2 + "/Clip_analysis/src/")
from metaexon_coverage import main_2h
sf_type_allowed = ['gc_rich_down', '', 'at_rich_down', 'ga_rich_down', 'ct_rich_down', 'other']
sys.path.insert(0, base1 + "/GC_AT_group_regulated_U1_U2/")
from spliceosome_regulation_enrichment import main_3g


@lp.parse(list_file="file", seddb="file", exon_type=["CCE", "ALL"],
          regulation=["up", "down"], output="dir", bed_folder="dir",
          refsize="file", metagene_script="file",
          sf_type=sf_type_allowed, reverse=["y", "n"])
def figure_creator(list_file, name_file, seddb, fasterdb, output, bed_folder,
                   refsize, metagene_script,
                   exon_type="CCE", regulation="down", nt="S",
                   sf_type=("", ""), reverse="n"):
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
    :param bed_folder: (str) a directory containing bed files
    :param refsize: (str) a file containing chromosome length of hg19 genome
    :param metagene_script: (str) path to metagene scriptS
    :param nt: (str) the nt we want to use for the figure 1.1D
    :param sf_type: (List(size=2, vtype=str)) the name of the 2 kind \
    of splicing factors groups of interest
    :param reverse: (str) reverse the figure 3g

    :return:
    """
    if output[-1] != "/":
        output += "/"
    print("Creating figure 1D...")
    main_1d(list_file.copy(), name_file.copy(), seddb, exon_type, regulation,
            output, nt)
    print("Creation of figure 2A...")
    main_2a(list_file.copy(), name_file.copy(), "C,S,A,T,G,W,Y,R".split(","),
            "2A_metaexon", exon_type, ["#0000FF", "#00aa00"], False,
            output, fasterdb)

    if sf_type[0] != "" and sf_type[1] != "":
        # Note, you have to modifiy the group factor.py script for other list
        fig_num = ["2.1B", "2.2B", "2.1C", "2.2C"]
        cols = ["iupac_downstream_intron_adjacent1"] * 2 + \
            ["iupac_upstream_intron_adjacent1"] * 2
        sf_type_tot = [sf_type, sf_type[::-1]] * 2
        ascending = [True, False] * 2
        for num, col, sf_type_cur, asc in zip(fig_num, cols, sf_type_tot,
                                              ascending):
            print("Creation of the figure %s" % num)
            main_2bc([col], "%s_%s_%s" % (num, col, sf_type_cur[0]),
                     seddb, exon_type, output, sf_type_cur[0], sf_type_cur[1],
                     regulation="down", contrast=20, operation="mean",
                     mascending=asc)
    print("Creating figure")
    main_2d(list_file.copy(), name_file.copy(), exon_type, output, seddb,
            fasterdb, fig_nums=("2.1D_", "2.2D_"))

    print("Creating figures 2E 2F and 2.1G")
    main_2efg(list_file.copy(), name_file.copy(),
              exon_type, seddb, fasterdb, output)

    print("Creation of the figure 2.2G")
    main_2g(list_file.copy(), name_file.copy(), exon_type, fasterdb,
            seddb, output)

    print("Creation of the figures 2H")
    main_2h(bed_folder, list_file.copy(), name_file.copy(), refsize, seddb,
            fasterdb, output, metagene_script, None, None, None, num_fig="2H")

    print("Creation of the figure 3G")
    dic_rev = {"y": True, "n": False}
    main_3g(list_file.copy(), name_file.copy(), seddb, output,
            dic_rev[reverse])

if __name__ == "__main__":
    figure_creator()