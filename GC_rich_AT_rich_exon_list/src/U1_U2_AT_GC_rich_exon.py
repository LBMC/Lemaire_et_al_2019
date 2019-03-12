#!/bin/usr/env python3


# -*- coding utf-8 -*-

"""
Description:

    The goal of this script is to test if GC_rich exon are more often regulated by U1 than U2. To do this we will \
    use a test frequency comparison.

    It also study the exons regulated by both factors regulating AT and GC rich exons. \
    To do this it display the fraction of common exon found in the list of every project done on factors \
    regulating AT and GC rich exons.

"""

# IMPORTS
import math
from ncephes import cprob
import union_dataset_function
import sqlite3
import os
import sys
import figure_creator
import numpy as np
import create_boxplot_cg_content
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from figure_creator import venn_diagram_creator

# FUNCTION
def exon_intersection(exon_list1, exon_list2):
    """
    Return the intersection of 2 exon lists

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the intersection between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon in exon_list2]
    return exon_list


def frequency_test(obs1, tot1, obs2, tot2):
    """
    Proportion test.

    :param obs1: (int) the count number of an amino acid X in the set of protein 1.
    :param tot1: (int) the total number of amino acids in the set of protein 1.
    :param obs2: (int) the count number of an amino acid X in the set of protein 2.
    :param tot2: (int) the total number of amino acids in the set of protein 2.
    :return: proportion test p-value
    """
    mean1 = float(obs1) / tot1
    mean2 = float(obs2) / tot2

    var1 = float(obs1) * (1 - mean1) * (1 - mean1) + (tot1 - obs1) * mean1 * mean1
    var2 = float(obs2) * (1 - mean2) * (1 - mean2) + (tot2 - obs2) * mean2 * mean2

    df = tot1 + tot2 - 2
    svar = (var1 + var2) / df
    t = (mean1-mean2) / math.sqrt(svar*(1.0/tot1 + 1.0/tot2))

    return cprob.incbet(0.5*df, 0.5, df/(df+t*t))


def get_all_exons_list(cnx, sf_list, sig, regulation):
    """
    Return every non-redundant exons regulated by at least one factor in ``sf_list`` (with the regulation \
    ``regulation``)

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param sig: (boolean) true to recover only significant splicing events (given by farline) or false \
    to recover every events
    :param regulation: (list of string) up or down or up + down, the regulation of event that we want \
    to recover
    :return: (list of list of int) list of exons shownig the regulation ``regulation`` at least for a factor \
    in ``sf_list``
    """
    exon_list = []
    for sf_name in sf_list:
        exon_list += union_dataset_function.get_events_4_a_sl_all(cnx, sf_name, regulation, sig)
    exon_list = union_dataset_function.washing_events_all(exon_list)
    return exon_list


def analysis_maker(gc_rich_exons, at_rich_exons, u1_exons, u2_exons, reg, sig, result_file):
    """
    Create a result file saying if GC_rich exons and AT rich exons are more often regulated by \
    U1 or U2.

    :param gc_rich_exons: (list of list of 2 int) list of GC rich exons identified by gene_id and exon_position on \
    gene
    :param at_rich_exons:(list of list of 2 int) list of AT rich exons identified by gene_id and exon_position on \
    gene
    :param u1_exons: (list of list of 2 int) list of exons regulated by U1
    :param u2_exons:(list of list of 2 int) list of exons regulated by U2
    :param reg: (string) up or down
    :param sig: (boolean) True if we want to recover every significant exon regulated, False to recover every exons \
    (significant and non-significant)
    :param result_file: (string) path where the result will be created
    """
    gc_rich_exon_u1_regulated = len(exon_intersection(gc_rich_exons, u1_exons))
    gc_rich_exon_u2_regulated = len(exon_intersection(gc_rich_exons, u2_exons))
    at_rich_exon_u1_regulated = len(exon_intersection(at_rich_exons, u1_exons))
    at_rich_exon_u2_regulated = len(exon_intersection(at_rich_exons, u2_exons))
    pvalue_gc = frequency_test(gc_rich_exon_u1_regulated, len(u1_exons),
                               gc_rich_exon_u2_regulated, len(u2_exons))
    pvalue_at = frequency_test(at_rich_exon_u1_regulated, len(u1_exons),
                               at_rich_exon_u2_regulated, len(u2_exons))
    if reg is None:
        reg = "all"
    if sig:
        name_sig = "significant"
    else:
        name_sig = "significant and non significant"
    with open(result_file, "a") as outfile:
        outfile.write("GC rich exons analysis with U1 and U2 %s %s-exons : U1 = %s / %s ; U2 = %s / %s - pvalue : %s\n"
                      % (name_sig, reg, gc_rich_exon_u1_regulated, len(u1_exons), gc_rich_exon_u2_regulated,
                         len(u2_exons), pvalue_gc))
        outfile.write("AT rich exons analysis with U1 and U2 %s %s-exons : U1 = %s / %s ; U2 = %s / %s  - pvalue : %s\n"
                      % (name_sig, reg, at_rich_exon_u1_regulated, len(u1_exons), at_rich_exon_u2_regulated,
                         len(u2_exons), pvalue_at))


def analysis_maker_projects(exon_list, name_project, u1_exons, u2_exons, reg, sig, result_file):
    """
    Create a result file saying if exons regulated in a project are more often regulated by \
    U1 or U2.

    :param exon_list: (list of list of 2 int) list of exon regulated in a particular project
    :param name_project:(string) list of AT rich exons identified by gene_id and exon_position on \
    gene
    :param u1_exons: (list of list of 2 int) list of exons regulated by U1
    :param u2_exons:(list of list of 2 int) list of exons regulated by U2
    :param reg: (string) up or down
    :param sig: (boolean) True if we want to recover every significant exon regulated, False to recover every exons \
    (significant and non-significant)
    :param result_file: (string) path where the result will be created
    """
    exon_u1_regulated = len(exon_intersection(exon_list, u1_exons))
    exon_u2_regulated = len(exon_intersection(exon_list, u2_exons))
    pvalue = frequency_test(exon_u1_regulated, len(u1_exons),
                            exon_u2_regulated, len(u2_exons))
    if pvalue < 0.05:
        if exon_u1_regulated / len(u1_exons) > exon_u2_regulated / len(u2_exons):
            enrichment = "U1 enriched"
        else:
            enrichment = "U2 enriched"
    else:
        enrichment = ""
    if reg is None:
        reg = "all"
    if sig:
        name_sig = "significant"
    else:
        name_sig = "significant and non significant"
    with open(result_file, "a") as outfile:
        outfile.write("%s exons analysis with U1 and U2 %s %s-exons : U1 = %s / %s ; U2 = %s / %s - pvalue : %s\t%s\n"
                      % (name_project, name_sig, reg, exon_u1_regulated, len(u1_exons), exon_u2_regulated,
                         len(u2_exons), pvalue, enrichment))


def extract_exon_list(filename):
    exon_list = []
    with open(filename, "r") as outfile:
        line = outfile.readline()
        while line:
            line = line.replace("\n", "")
            line = line.split("\t")
            exon_list.append(line)
            line = outfile.readline()
    return exon_list


def difference(exon_list1, exon_list2):
    """
    Return the difference of 2 exon lists

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the intersection between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon not in exon_list2]
    return exon_list


def get_delta_psi(cnx, exon, id_projects, sig):
    """
    Get the delta psi of exons regulated by this category of splicing factor

    :param cnx: (sqlite3 connect object) connection to seddb
    :param exon: (list of strings and 2 int) list of exons
    :param id_projects: (tuple of int) tuple of id projects
    :param sig: (boolean) True if you want only significant exons False else
    :return: (list of list of one string and 1 int) list of exons with the deltapsi associated if it's regulated by \
     ``category`` of factor
    """
    cursor = cnx.cursor()
    delta_values = []
    for my_id in id_projects:
        if not sig:
            query = """
                    SELECT DISTINCT t2.project_name, t1.delta_psi
                    FROM ase_event t1, rnaseq_projects t2
                    WHERE t2.id = t1.id_project
                    AND t2.id = %s
                    AND t1.gene_id = %s
                    AND t1.exon_skipped = %s
                    AND t1.delta_psi IS NOT NULL""" % (my_id, exon[0], exon[1])
            cursor.execute(query)
            res = cursor.fetchall()
        else:
            cursor = cnx.cursor()
            query = """
                    SELECT DISTINCT t2.project_name, t1.delta_psi
                    FROM ase_event t1, rnaseq_projects t2
                    WHERE t2.id = t1.id_project
                    AND t2.id = %s
                    AND t1.gene_id = %s
                    AND t1.exon_skipped = %s
                    AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
                    AND pvalue_glm_cor <= 0.05""" % (my_id, exon[0], exon[1])
            cursor.execute(query)
            res = cursor.fetchall()
            if len(res) == 0:
                query = """
                        SELECT DISTINCT t2.project_name, t1.delta_psi
                        FROM ase_event t1, rnaseq_projects t2
                        WHERE t2.id = t1.id_project
                        AND t2.id = %s
                        AND t1.gene_id = %s
                        AND t1.exon_skipped = %s
                        AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
                        AND pvalue <= 0.05""" % (my_id, exon[0], exon[1])
                cursor.execute(query)
                res = cursor.fetchall()
        if res:
            if len(res) > 1:
                print("WARNING : 2 detapsi where retrieved yet one project was given")
            delta_values.append(res[0])
    if delta_values:
        values = []
        for val in delta_values:
            values.append(val[1])
        mean_delta_psi = np.mean(values)
        res_exon = [["%s_%s" % (exon[0], exon[1]), mean_delta_psi]]
        return res_exon
    return [["%s_%s" % (exon[0], exon[1]), None]]


def get_list_project(cnx, sf_list):
    """
    get the list of project id that was done on the ``category`` of factor chosen

    :param cnx: (sqlite3 connect object) connection to seddb
    :param sf_list: (strings) list of interest sf
    :return: (list of int) list of id project
    """
    cursor = cnx.cursor()
    id_projects = []
    for sf_name in sf_list:
        query = """SELECT DISTINCT id from rnaseq_projects WHERE sf_name = '%s'""" % sf_name
        cursor.execute(query)
        res = cursor.fetchall()
        for my_id in res:
            id_projects.append(my_id[0])
    return tuple(id_projects)


def get_the_deltapsi_associated_to_a_categorie_of_factor(cnx, exon_list, sf_control, sig):
    """
    Get the deltapsi associated to a list of splicing factor in a particular cell line regulated by a particular \
    category of splicing factor.

    :param cnx: (sqlite3 connect object) connection to seddb
    :param exon_list: (list of list of strings and 2 int) list of exons
    :param sf_control: (string) list of control sf
    :param sig: (boolean) True if you want to get only significant exon False else
    :return: (dictionary) a dictionary of 2 keys, the first keys links to a list containing the exon names, \
    the second key links to a list where associated pvalue are stored
    """
    id_projects = get_list_project(cnx, sf_control)
    # new_t = []
    # for val in exon_list_control:
    #     new_t.append("%s_%s" % (get_gene_name(cnx, val[0]) , val[1]))
    # new_t = sorted(new_t)
    # for val in new_t:
    #     print(val)
    # exon_list = exon_list[0:100]
    exon_dic_completed = {"name": [], "delta_psi": []}
    count = 0
    for exon in exon_list:
        count += 1
        sys.stdout.write("%s / %s \r" % (count, len(exon_list)))
        sys.stdout.flush()
        exon_completed = get_delta_psi(cnx, exon, id_projects, sig)
        if len(exon_completed) == 1:
            if exon_completed[0][-1] is not None:
                exon_dic_completed["name"].append(exon_completed[0][0])
                exon_dic_completed["delta_psi"].append(exon_completed[0][1])
        else:
            for i in range(len(exon_completed)):
                exon_dic_completed["name"].append(exon_completed[i][0])
                exon_dic_completed["delta_psi"].append(exon_completed[i][1])
    return exon_dic_completed


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def create_figure(list_values, list_name, output, name_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    color_list = ['#1f77b4', '#2ca02c', '#1f77b4', '#2ca02c', '#FF8000',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    pval_u1 = mann_withney_test_r(list_values[0], list_values[1])
    pval_u2 = mann_withney_test_r(list_values[2], list_values[3])
    pval_gc = mann_withney_test_r(list_values[0], list_values[2])
    pval_at = mann_withney_test_r(list_values[1], list_values[3])
    data = []
    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, # "fillcolor": color_list[i], "opacity": 0.6,
                     "line": {"color": color_list[i]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

    layout = go.Layout(
        title='deltapsi of GC/AT exons in U1 or U2 experiments '
              '<br> mann whitney test U1 AT exons vs GC exons : p = %.2E'
              '<br> mann whitney test U2 AT exons vs GC exons : p = %.2E'
              '<br> mann whitney test GC U1 exons vs U2 exons : p = %.2E'
              '<br> mann whitney test AT U1 exons vs U2 exons : p = %.2E'
              % (pval_u1, pval_u2, pval_gc, pval_at),
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
    plotly.offline.plot(fig, filename="%s%s_deltapsi.html" % (output, name_fig),
                    auto_open=False, validate=False)


def main():
    """
    Test if GC_rich exon are more often regulated by U1 than U2. To do this we will \
    use a test frequency comparison. \
    It also studies the exons regulated by both factors regulating AT and GC rich exons. \
    To do this it display the fraction of common exon found in the list of every project done on factors \
    regulating AT and GC rich exons.

    """
    debug = 0  # disabled
    regulation = "down"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/")
    print("test")
    result_file = output + "result.txt"
    at_file = "%sAT_rich_exons" % output
    gc_file = "%sGC_rich_exons" % output
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    list_file = [gc_file, at_file]
    exon_list_files = []
    for my_file in list_file:
        exon_list_files.append(extract_exon_list(my_file))

    gc_exon_list = exon_list_files[0]
    at_exon_list = exon_list_files[1]
    # write_origin(intersect, gc_and_at_rich_exon, cnx, regulation)
    with open(result_file, "w") as outfile:
        outfile.write("GC_rich_exon : %s - AT_rich exons : %s\n" %
                      (len(exon_list_files[0]), len(exon_list_files[1])))


    # u1_exons = get_all_exons_list(cnx, figure_creator.u1_factors, True, regulation)
    # u2_exons = get_all_exons_list(cnx, figure_creator.u2_factors, True, regulation)
    # common_u1_u2 = exon_intersection(u1_exons, u2_exons)
    # u1_pur = difference(u1_exons, common_u1_u2)
    # print("U1 sig : %s " % len(u1_pur))
    # u1_exons_all = get_all_exons_list(cnx, figure_creator.u1_factors, False, regulation)
    # # u2_exons_all = get_all_exons_list(cnx, figure_creator.u2_factors, False, regulation)
    # # common_u1_u2 = exon_intersection(u1_exons_all, u2_exons_all)
    # # u1_pur_all = difference(u1_exons_all, common_u1_u2)
    # print("U1 non sig : %s " % len(u1_exons_all))
    # print("interestction u1 sig non sig : %s" % len(exon_intersection(u1_pur, u1_exons_all)))
    significant = [True, False]
    regulations = ["down", "all"]
    for regulation in regulations:
        for sig in significant:
            print(str(sig), str(regulation))
            u1_exons = get_all_exons_list(cnx, figure_creator.u1_factors, sig, regulation)
            print(len(u1_exons))
            u2_exons = get_all_exons_list(cnx, figure_creator.u2_factors, sig, regulation)
            print(len(u2_exons))

            common_u1_u2 = exon_intersection(u1_exons, u2_exons)
            u1_pur = difference(u1_exons, common_u1_u2)
            u2_pur = difference(u2_exons, common_u1_u2)
            analysis_maker(gc_exon_list, at_exon_list, u1_pur, u2_pur, regulation, sig, result_file)
            if sig:
                name_sig = "sig"
            else:
                name_sig = "non_sig"
            venn_diagram_creator(u1_exons, "U1_%s_%s" % (name_sig, regulation), u2_exons, "U2_%s_%s" % (name_sig, regulation), output)

            exon_u1_gc = exon_intersection(gc_exon_list, u1_pur)
            print("exon_u1_gc %s" % len(exon_u1_gc))
            u1_gc_dic = get_the_deltapsi_associated_to_a_categorie_of_factor(cnx, exon_u1_gc, figure_creator.u1_factors, sig)
            exon_u1_at = exon_intersection(at_exon_list, u1_pur)
            print("exon_u1_at %s" % len(exon_u1_at))
            u1_at_dic = get_the_deltapsi_associated_to_a_categorie_of_factor(cnx, exon_u1_at, figure_creator.u1_factors, sig)

            exon_u2_gc = exon_intersection(gc_exon_list, u2_pur)
            print("exon_u2_gc %s" % len(exon_u2_gc))
            u2_gc_dic = get_the_deltapsi_associated_to_a_categorie_of_factor(cnx, exon_u2_gc, figure_creator.u2_factors, sig)
            exon_u2_at = exon_intersection(at_exon_list, u2_pur)
            print("exon_u2_at %s" % len(exon_u2_at))
            u2_at_dic = get_the_deltapsi_associated_to_a_categorie_of_factor(cnx, exon_u2_at, figure_creator.u2_factors, sig)

            name_group = ["u1_gc_deltapsi", "u1_at_deltapsi", "u2_gc_deltapsi", "u2_at_deltapsi"]
            list_val = [u1_gc_dic["delta_psi"], u1_at_dic["delta_psi"], u2_gc_dic["delta_psi"], u2_at_dic["delta_psi"]]
            create_figure(list_val, name_group, output, "%s_%s" % (name_sig, regulation))


if __name__ == "__main__":
    main()
