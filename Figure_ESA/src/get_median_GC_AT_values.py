#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    Get the median GC content of every exons regulated by each splicing factor of interest
"""

#IMPORTS
import correlation_maker
import os
import sqlite3
import group_factor
import union_dataset_function
import numpy as np
import pandas as pd
import figure_producer
import sys


def main(level):
    """
    Create the correlation matrix (gene_size vs iupac)
    """
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = sqlite3.connect(seddb)
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/gc_at_median/"
    dic_res = {"name": [], "up": [], "down": []}
    if level == "project":
        name_fig = "median_GC_AT_project.txt"
        id_projects, name_projects = group_factor.get_id_and_name_project_wanted(cnx, None)
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for regulation in ["up", "down"]:
            for i in range(len(id_projects)):
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
                dic_res["name"].append(name_projects[i])
                values = correlation_maker.get_interest_values(cnx, exon_list, "iupac_exon", "S")
                dic_res[regulation].append(np.nanmedian(values))
    else:
        name_fig = "median_GC_AT_sf.txt"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        sf_list = group_factor.get_wanted_sf_name(None)
        for sf_name in sf_list:
            dic_res["name"].append(sf_name)
            for regulation in ["up", "down"]:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
                values = correlation_maker.get_interest_values(cnx, exon_list, "iupac_exon", "S")
                dic_res[regulation].append(round(np.nanmedian(values), 1))
    df = pd.DataFrame(dic_res)
    df = df[["name", "up", "down"]]
    df.to_csv(output + name_fig, sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1])