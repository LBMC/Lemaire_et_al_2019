#!/usr/bin/python3

# -*- coding utf-8 -*-


import os
import sqlite3
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import numpy as np
import math
import sys
import group_factor
import union_dataset_function


def get_control_GC_content(cnx, exon_type):
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT iupac_exon
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT iupac_exon
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    gc_content = []
    for iupac in tuple_list:
        gc_content.append(iupac[0].split(";")[4])
    # turn tuple into list
    return gc_content


def get_control_median_flanking_intron_size(cnx, exon_type):
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT upstream_intron_size, downstream_intron_size
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT upstream_intron_size, downstream_intron_size
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    median_intron_size = []
    for size in tuple_list:
        median_intron_size.append(np.nanmedian(np.array([size[0], size[1]], dtype=float)))
    # turn tuple into list
    return median_intron_size


def calculate_gc_content(cnx, gene_id, exon_pos):
    cursor = cnx.cursor()
    query = "SELECT iupac_exon FROM sed WHERE gene_id = ? and exon_pos = ?"
    cursor.execute(query, (gene_id, exon_pos,))
    res = cursor.fetchone()
    return res[0].split(";")[4]


def calculate_median_flanking_intron_size(cnx, gene_id, exon_pos):
    cursor = cnx.cursor()
    query = "SELECT upstream_intron_size, downstream_intron_size FROM sed WHERE gene_id = ? and exon_pos = ?"
    cursor.execute(query, (gene_id, exon_pos,))
    res = cursor.fetchone()
    return np.nanmedian(np.array([res[0], res[1]], dtype=float))


def extract_gc_content_from_file(cnx, filename):
    list_gc = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            list_gc.append(calculate_gc_content(cnx, line[0], line[1]))
            line = in_file.readline()
    return list_gc

def extract_gc_content_from_list(cnx, exon_list):
    list_gc = []
    for exon in exon_list:
        list_gc.append(calculate_gc_content(cnx, exon[0], exon[1]))
    return list_gc

def extract_median_flanking_intron_size_from_file(cnx, filename):
    median_intron_size = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            median_intron_size.append(calculate_median_flanking_intron_size(cnx, line[0], line[1]))
            line = in_file.readline()
    return median_intron_size



def create_figure(list_values, list_name, output, regulation, name_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    color_list=['#1f77b4', '#2ca02c', '#1f77b4', '#2ca02c', '#9467bd',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    data = []
    if "U1" in name_fig:
        color_list = ['navy', 'green', 'purple', 'red', 'navy',
                      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        pval1 = mann_withney_test_r(list_values[0], list_values[3])
        pval2 = mann_withney_test_r(list_values[1], list_values[3])
        pval3 = mann_withney_test_r(list_values[2], list_values[3])
        title = """%s of %s exons regulated by different factors
              <br> mw test %s exons vs CCE exons : p = %.2E, mw test %s exons vs CCE exons : p = %.2E
              <br> mw test %s exons vs CCE exons : p = %.2E""" % (name_fig, regulation, list_name[0], pval1, list_name[1], pval2, list_name[2], pval3)
    if "U2" in name_fig:
        color_list = ['navy', 'green', 'purple', '#e377c2', 'navy',
                      '#8c564b', 'red', '#7f7f7f', '#bcbd22', '#17becf']
        values = []
        for i in range(len(list_values) - 1):
            values.append(list_name[i])
            values.append(list_name[-1])
            values.append(mann_withney_test_r_lesser(list_values[i], list_values[-1]))
        title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)
        tmp = ""
        for i in range(0, len(values), 3):
            tmp += "<br> %s vs %s, p=%.2E" % (values[i], values[i + 1], values[i + 2])
        title += tmp
    elif "GC" in name_fig:
        pval_at_cg = mann_withney_test_r(list_values[0], list_values[1])
        pval_u1_u2 = mann_withney_test_r(list_values[2], list_values[3])
        title = """%s of %s exons regulated by different factors
        <br> mann whitney test AT exons vs GC exons : p = %.2E
        <br> mann whitney test U1 exons vs U2 exons : p = %.2E""" % (name_fig, regulation, pval_at_cg, pval_u1_u2)
    elif  "median flanking intron size" in name_fig:
        name_fig = "log " + name_fig
        for i in range(len(list_values)):
            list_values[i] = list(map(math.log, list_values[i]))

        pval_at_cg = mann_withney_test_r_lesser(list_values[0], list_values[1])
        pval_u1_u2 = mann_withney_test_r_lesser(list_values[2], list_values[3])
        title = """%s of %s exons regulated by different factors
                <br> mann whitney test AT exons vs GC exons : p = %.2E
                <br> mann whitney test U1 exons vs U2 exons : p = %.2E""" % (
        name_fig, regulation, pval_at_cg, pval_u1_u2)
    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "fillcolor": color_list[i], "opacity": 0.6,
                     "line": {"color": "black"},
                     "box": {"visible": True}, "meanline": {"visible": True}})

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title=name_fig,
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
    plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl.html" % (output, name_fig, regulation),
                    auto_open=False, validate=False)



def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='greater', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval

def mann_withney_test_r_lesser(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='less', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def main():
    exon_type = "CCE"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    if len(sys.argv) < 2 or sys.argv[1] == "pure":
        u1_file = "%sU1_exons" % output
        u2_file = "%sU2_exons" % output
        at_file = "%sAT_rich_exons" % output
        gc_file = "%sGC_rich_exons" % output
    else:
        u1_file = "%sU1_with_intersection_exons" % output
        u2_file = "%sU2_with_intersection_exons" % output
        at_file = "%sAT_rich_with_intersection_exons" % output
        gc_file = "%sGC_rich_with_intersection_exons" % output
    name_file = ["GC_rich_exons", "AT_rich_exons", "u1_exons", "u2_exons", exon_type]
    list_file = [gc_file, at_file, u1_file, u2_file, None]
    list_gc_content = []
    list_median_flanking_intron_size = []
    for i in range(len(name_file)):
        if name_file[i] != exon_type:
            list_gc_content.append(extract_gc_content_from_file(cnx, list_file[i]))
        else:
            list_gc_content.append(get_control_GC_content(cnx, exon_type))
    create_figure(list_gc_content, name_file, output, "down", "GC content " + sys.argv[1])

    for i in range(len(name_file)):
        if name_file[i] != exon_type:
            list_median_flanking_intron_size.append(extract_median_flanking_intron_size_from_file(cnx, list_file[i]))
        else:
            list_median_flanking_intron_size.append(get_control_median_flanking_intron_size(cnx, exon_type))
    create_figure(list_median_flanking_intron_size, name_file, output, "down", "median flanking intron size " + sys.argv[1])
    # --  same thing but with U1 factors -----
    sf_list = list(group_factor.u1_factors) + ['DDX5_DDX17'] + ["CCE"]
    list_gc_content = []
    for sf_name in sf_list:
        if sf_name != exon_type:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, "down")
            list_gc_content.append(extract_gc_content_from_list(cnx, exon_list))
            print("%s : %s exons" % (sf_name, len(exon_list)))
        else:
            list_gc_content.append(get_control_GC_content(cnx, exon_type))
    create_figure(list_gc_content, sf_list, output, "down", "GC content U1")
    # --  same thing but with U2 factors -----
    sf_list = list(group_factor.u2_factors)  + ["CCE"]
    list_gc_content = []
    for sf_name in sf_list:
        if sf_name != exon_type:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, "down")
            list_gc_content.append(extract_gc_content_from_list(cnx, exon_list))
            print("%s : %s exons" % (sf_name, len(exon_list)))
        else:
            list_gc_content.append(get_control_GC_content(cnx, exon_type))
    create_figure(list_gc_content, sf_list, output, "down", "GC content U2")


if __name__ == "__main__":
    main()