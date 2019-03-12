#!/usr/bin/python3

# -*- coding: utf-8 -*-


# from figure_creator import at_rich_down, gc_rich_down, u1_factors, u2_factors
import pymysql
import conf
import sqlite3
import os
import numpy as np
import sys
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v

## various function useful for any other script
## ( made to avoid circular imports )


def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = pymysql.connect(user=conf.user, password=conf.password, host=conf.host, database=conf.database)
    return cnx



def get_psi_control_of_an_exon(cnx, gene_symbol, coordinate, exon_pos):
    cursor = cnx.cursor()
    query = """SELECT PSI FROM PSI_lignees_Controle
           WHERE gene_symbole = '%s'
           AND exon = %s
           AND coordonnees = '%s'""" % (gene_symbol, exon_pos, coordinate)
    cursor.execute(query)
    res = cursor.fetchall()
    list_val = [value[0] for value in res]
    if len(list_val) == 0:
        return None
    return np.nanmedian(np.array(list_val, dtype=float))


def get_dpsi_control_of_an_exon(cnx, gene_symbol, coordinate, exon_pos):
    cursor = cnx.cursor()
    query = """SELECT PSI FROM PSI_lignees_Controle
           WHERE gene_symbole = '%s'
           AND exon = %s
           AND coordonnees = '%s'""" % (gene_symbol, exon_pos, coordinate)
    cursor.execute(query)
    res = cursor.fetchall()
    list_val = [value[0] for value in res]
    list_val = np.array(list_val, dtype=float)
    list_val = list_val[~np.isnan(list_val)]
    if len(list_val) < 2:
        return None
    median = np.nanmedian(list_val)
    var = 0
    for val in list_val:
       var +=  (val - median) ** 2
    variance = var / len(list_val)
    return variance


def get_gene_symbol_exon(cnx, gene_id, exon_pos):
    cursor = cnx.cursor()
    query = """SELECT DISTINCT gene_symbol, chromosome, start, stop FROM ase_event
               WHERE gene_id = %s
               AND exon_skipped = %s""" % (gene_id, exon_pos)
    cursor.execute(query)
    res = cursor.fetchall()
    if len(res) > 1:
        print("More than one gene_symbol retrieved for an exon")
        exit(1)
    return res[0][0], "%s:%s-%s" % (res[0][1], res[0][2], res[0][3])


def extract_median_delta_psi_from_file(cnx, sf_cnx, filename):
    list_psi = []
    with open(filename, "r") as in_file:
        print("Collecting psi from %s file " % filename)
        count = 0
        line = in_file.readline()
        while line:
            count += 1
            sys.stdout.write("%s    \r" % count)
            sys.stdout.flush()
            line = line.replace("\n", "")
            line = line.split("\t")
            gene_symbol, coordinates = get_gene_symbol_exon(cnx, line[0], line[1])
            list_psi.append(get_dpsi_control_of_an_exon(sf_cnx, gene_symbol, coordinates, line[1]))
            line = in_file.readline()
    list_psi = np.array(list_psi, dtype="float")
    return list(list_psi[~np.isnan(list_psi)])


def extract_median_psi_from_file(cnx, sf_cnx, filename):
    list_psi = []
    with open(filename, "r") as in_file:
        print("Collecting psi from %s file " % filename)
        count = 0
        line = in_file.readline()
        while line:
            count += 1
            sys.stdout.write("%s    \r" % count)
            sys.stdout.flush()
            line = line.replace("\n", "")
            line = line.split("\t")
            gene_symbol, coordinates = get_gene_symbol_exon(cnx, line[0], line[1])
            list_psi.append(get_psi_control_of_an_exon(sf_cnx, gene_symbol, coordinates, line[1]))
            line = in_file.readline()
    list_psi = np.array(list_psi, dtype="float")
    return list(list_psi[~np.isnan(list_psi)])


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def create_figure(list_values, list_name, output, regulation, name_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    data = []
    color_list=['#1f77b4', '#2ca02c', '#1f77b4', '#2ca02c', '#9467bd',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    pval_at_cg = mann_withney_test_r(list_values[0], list_values[1])
    pval_u1_u2 = mann_withney_test_r(list_values[2], list_values[3])
    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "line": {"color": color_list[i]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

    layout = go.Layout(
        title='%s of %s exons regulated by different factors '
              '<br> mann whitney test AT exons vs GC exons (two sided) : p = %.2E'
              '<br> mann whitney test U1 exons vs U2 exons (two sided): p = %.2E' % (name_fig, regulation, pval_at_cg, pval_u1_u2),
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

def main():
    exon_type = "CCE"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    sl_cnx = connection()
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
    name_file = ["GC_rich_exons", "AT_rich_exons", "u1_exons", "u2_exons"]
    list_file = [gc_file, at_file, u1_file, u2_file]
    list_psi = []
    list_dpsi = []
    for i in range(len(name_file)):
        list_psi.append(extract_median_psi_from_file(cnx, sl_cnx,list_file[i]))
        print(list_psi)
    create_figure(list_psi, name_file, output, "down", "PSI " + sys.argv[1])

    for i in range(len(name_file)):
       list_dpsi.append(extract_median_delta_psi_from_file(cnx, sl_cnx, list_file[i]))
    create_figure(list_dpsi, name_file, output, "down", "var psi" + sys.argv[1])


if __name__ == "__main__":
    main()