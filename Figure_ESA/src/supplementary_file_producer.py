#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    This file allows to create 3 figures displaying for each splicing factor their GC content in \
    exon and their surrounding intronic sequence of 100 nt.
"""

import figure_producer
import os
import sys
import group_factor
import numpy as np
import plotly.graph_objs as go
import plotly
import math

nt_list = ["S"]


def create_figure_iupac_dnt(cnx, name_projects, target_column, regulation, output, nt_dnt):
    """
    Create a figure for every column in sed database whose name contains "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param name_projects: (list of string) the list of sf_name
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    :param nt_dnt: (string) the nt of interest or the di-nucleotide of interest
    """
    if "iupac" in target_column:
        result = figure_producer.get_values_for_many_projects_iupac_dnt(cnx, name_projects, target_column, regulation, nt_dnt, True)
    else:
        result = figure_producer.get_values_for_many_projects(cnx, name_projects, target_column, regulation, True)

    target_column = target_column.replace("iupac", "%s_nt" % nt_dnt)
    d = {name_projects[i]: result[i] for i in range(len(result))}
    e = sorted(d.items(), key=lambda x: np.median(x[1]), reverse=True)
    new_name = [x[0] for x in e]
    new_result = [x[1] for x in e]
    log = False
    for list_res in new_result:
        if max(list_res) > 1000:
            log = True
    if log:
        new_result = [list(map(math.log10, x[1])) for x in e]

    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(new_result))]
    data = []
    for i in range(len(new_result)):
        data.append({"y": new_result[i], "type": "violin",
                     "name": new_name[i], "visible": True, "marker": {"color": c[i]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

    layout = go.Layout(
        title='%s of %s exons for every projects' % (target_column, regulation),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_exons_figure.html" % (output, target_column, regulation),
                        auto_open=False, validate=False)


def main():
    """
    Launch the creation of figures.
    """
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    regs = ["down"]
    cnx = figure_producer.connexion(seddb)
    columns = ["iupac_exon", "iupac_upstream_intron_proxi", "iupac_downstream_intron_proxi", "median_flanking_intron_size"]
    if len(sys.argv) < 2:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/supplementary_figure/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        sf_list = group_factor.get_wanted_sf_name(None)
        for regulation in regs:
            print(regulation)
            for target_column in columns:
                print("   %s" % target_column)
                if "iupac" in target_column:
                    for nt in nt_list:
                        create_figure_iupac_dnt(cnx, sf_list, target_column, regulation, output, nt)
                else:
                    create_figure_iupac_dnt(cnx, sf_list, target_column, regulation, output, None)


if __name__ == "__main__":
    main()
