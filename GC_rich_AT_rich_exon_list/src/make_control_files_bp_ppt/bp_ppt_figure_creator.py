#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    This script will create 2 figures one displaying the ppt score of 2 the GC and at list of exons
"""

import os
import sqlite3
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import sys
import exon_class
import function
sys.path.insert(0, "/media/nicolas/DD_1/Splicing_Lore_project/GC_rich_AT_rich_exon_list/src")
import group_factor
import union_dataset_function
import numpy as np


def get_redundant_list_of_value(cnx, exon_list, target_column, with_none=False):
    """
    Get the individual values for ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of ExonClass object)list of exons
    :param target_column: (string) the column for which we want to get information on exons.
    :param with_none: (boolean) False is we want no none True else
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    for exon in exon_list:
        query = """SELECT %s
                   FROM sed
                   where gene_id = %s
                   AND exon_pos = %s """ % (target_column, exon.gene.id, exon.position)
        cursor.execute(query)
        r = cursor.fetchone()[0]
        if r is not None:
            res.append(r)
        elif r is None and with_none:
            res.append(None)
    return res


def extract_exon_files(cnx, filename):
    """
    :param cnx: (sqlite3 connect object) connection to fasterDB lite
    :param filename: (string) the name of a file containing exons
    :return: (list of Exonclass object) list of exons
    """
    exon_list = []
    with open(filename, "r") as outfile:
        line = outfile.readline()
        while line:
            line = line.replace("\n", "")
            line = line.split("\t")
            exon = exon_class.ExonClass(cnx, str(line[0]), int(line[0]), int(line[1]))
            exon_list.append(exon)
            line = outfile.readline()
    return exon_list


def get_control_exon_information(cnx, exon_type, target_column):
    """
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param target_column:  (string) the column of interest
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT {}
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(target_column, exon_type)
    else:
        query = """SELECT {}
                   FROM sed
                """.format(target_column)
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    values = []
    for val in result:
        values.append(val[0])
    return values


def create_figure(list_values, list_name, output, regulation, name_fig, type_fig="group"):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of figure group or U1
    """
    color_list=['#1f77b4', '#1fa7b4', '#2ca02c',  '#2ca05c', 'red',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    data = []
    if type_fig == "group":
        pval_gc_pure = mann_withney_test_r(list_values[0], list_values[4])
        pval_gc_all = mann_withney_test_r(list_values[1], list_values[4])
        pval_at_pure = mann_withney_test_r(list_values[2], list_values[4])
        pval_at_all = mann_withney_test_r(list_values[3], list_values[4])
        title = """%s of %s exons regulated by different factors
              <br> mw test GC pure vs CCE, p=%.2E, GC all vs CCE, p = %.2E
              <br> mw test AT pure vs CCE, p=%.2E, AT all vs CCE, p = %.2E""" % (name_fig, regulation, pval_gc_pure, pval_gc_all, pval_at_pure, pval_at_all)
    else:
        color_list = ['navy', 'green', 'purple', '#e377c2', 'navy',
                      '#8c564b', 'red', '#7f7f7f', '#bcbd22', '#17becf']
        values = []
        for i in range(len(list_values) - 1):
            values.append(list_name[i])
            values.append(list_name[-1])
            values.append(mann_withney_test_r(list_values[i], list_values[-1]))
        title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)
        tmp = ""
        for i in range(0, len(values), 3):
            tmp += "<br> %s vs %s, p=%.2E" % (values[i], values[i + 1], values[i + 2])
        title += tmp

    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, # "fillcolor": color_list[i], "opacity": 0.6,
                     "line": {"color": color_list[i]},
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


def create_barplot(list_values, list_name, output, regulation, name_fig):
    """
    Create a barplot showing the values of ``name_fig`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    if len(list_values) < 6:
        color_list=['#1f77b4', '#1fa7b4', '#2ca02c',  '#2ca05c', 'red']
    else:
        color_list=['#03b4d9','#0000FF', 'purple', 'green', 'olive', 'navy', 'pink', 'orange', "red"]
    data = []
    abscissa = ["<=1", "=2", ">=3"]
    for i in range(len(list_values)):
        ordinate = []
        ordinate.append(sum(np.array(list_values[i]) <= 1) / len(list_values[i]))
        ordinate.append(sum(np.array(list_values[i]) == 2) / len(list_values[i]))
        ordinate.append(sum(np.array(list_values[i]) >= 3) / len(list_values[i]))
        data.append(go.Bar(x=abscissa, y=ordinate, name=list_name[i], marker=dict(color=color_list[i])))

    title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)

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
        xaxis=dict(title="nb branch point"),
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
    plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl_barplot.html" % (output, name_fig, regulation),
                    auto_open=False, validate=False)


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def get_bp_ppt_score_list(output, exon_list, name_list, size):
    """
    Calculate or retrieve bp and ppt score
    :param output: (string) path where the result will be stored or retrieved
    :param exon_list: (list of ExonClass object) List of exons
    :param name_list: (string) the name of the exons list
    :param size: (int) the size of window wanted
    :return: (2 list of floats) list of bp score and list of ppt score
    """
    name_store_file = "%s%s_%s_bp_ppt_score.py" % (output, name_list, size)
    if not os.path.isfile(name_store_file):
        print("Calculating ppt bp score using SVM BP FINDER")
        bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list = function.bp_ppt_calculator(exon_list, size)
        with open(name_store_file, "w") as outfile:
            outfile.write("bp_score_list=%s\n" % str(bp_score_list))
            outfile.write("ppt_score_list=%s\n" % str(ppt_score_list))
            outfile.write("nb_bp_list=%s\n" % str(nb_bp_list))
            outfile.write("nb_good_bp_list=%s" % str(nb_good_bp_list))
    else:
        print("recovering ppt and bp score already stored in %s" % name_store_file)
        sys.path.insert(0, output)
        stored = __import__("%s_%s_bp_ppt_score" % (name_list, size))
        bp_score_list = stored.bp_score_list
        ppt_score_list = stored.ppt_score_list
        nb_bp_list = stored.nb_bp_list
        nb_good_bp_list = stored.nb_good_bp_list
    return bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list


def main():
    exon_class.set_debug(0)
    exon_type = "CCE"
    ctrl_output = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt", "result/make_control_files_bp_ppt/")
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt",
                                                                      "result/bp_ppt_score/")
    file_dir = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt", "result/")
    if not os.path.isdir(ctrl_output):
        os.mkdir(ctrl_output)
    if not os.path.isdir(output):
        os.mkdir(output)
    ctrl_dir = os.path.realpath(os.path.dirname(__file__)) + "/control_dictionaries/"
    sys.path.insert(0, ctrl_dir)
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt", "data/fasterDB_lite.db")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt", "data/sed.db")
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)

    at_pure_file = "%sAT_rich_exons" % file_dir
    gc_pure_file = "%sGC_rich_exons" % file_dir
    at_all_file = "%sAT_rich_with_intersection_exons" % file_dir
    gc_all_file = "%sGC_rich_with_intersection_exons" % file_dir
    name_file = ["GC_pure_exons", "GC_all_exons", "AT_pure_exons", "AT_all_exons", exon_type]
    list_file = [gc_pure_file, gc_all_file, at_pure_file, at_all_file, None]

    dict_score_3ss = {i:{"bp_score_list":[], "ppt_score_list":[], "nb_bp_list":[], "nb_good_bp_list":[]} for i in [100, 50, 25]}
    list_force_acceptor = []
    list_force_donor = []
    for i in range(len(name_file)):
        if name_file[i] != exon_type:
            exon_list = extract_exon_files(cnx, list_file[i])
            list_force_acceptor.append(get_redundant_list_of_value(cnx_sed, exon_list, "force_acceptor"))
            list_force_donor.append(get_redundant_list_of_value(cnx_sed, exon_list, "force_donor"))
            for size in dict_score_3ss.keys():
                bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list = get_bp_ppt_score_list(ctrl_output, exon_list, name_file[i], size)
                dict_score_3ss[size]["bp_score_list"].append(bp_score_list)
                dict_score_3ss[size]["ppt_score_list"].append(ppt_score_list)
                dict_score_3ss[size]["nb_bp_list"].append(nb_bp_list)
                dict_score_3ss[size]["nb_good_bp_list"].append(nb_good_bp_list)

        else:
            for size in dict_score_3ss.keys():
                mod = __import__("%s_%s_bp_ppt_score" % (exon_type, size))
                dict_score_3ss[size]["bp_score_list"].append(mod.bp_score)
                dict_score_3ss[size]["ppt_score_list"].append(mod.ppt_score)
                dict_score_3ss[size]["nb_bp_list"].append(mod.nb_bp)
                dict_score_3ss[size]["nb_good_bp_list"].append(mod.nb_good_bp)
            list_force_acceptor.append(get_control_exon_information(cnx_sed, exon_type, "force_acceptor"))
            list_force_donor.append(get_control_exon_information(cnx_sed, exon_type, "force_donor"))

    for size in dict_score_3ss.keys():
        print("------------> %s nt " % size)
        # create_figure(dict_score_3ss[size]["ppt_score_list"], name_file, output, "down", "ppt_score_(%snt)" % size)
        # create_figure(dict_score_3ss[size]["bp_score_list"], name_file, output, "down", "bp_score_(%snt)" % size )
        # create_figure(dict_score_3ss[size]["nb_bp_list"], name_file, output, "down", "nb_branch_point_(%snt)" % size)
        # create_figure(dict_score_3ss[size]["nb_good_bp_list"], name_file, output, "down", "nb_good_branch_point_(%snt)" % size)
        create_barplot(dict_score_3ss[size]["nb_bp_list"], name_file, output, "down", "prop_nb_branch_point_(%snt)" % size)
        create_barplot(dict_score_3ss[size]["nb_good_bp_list"], name_file, output, "down", "prop_nb_good_branch_point_(%snt)" % size)
    create_figure(list_force_acceptor, name_file, output, "down", "force_acceptor")
    create_figure(list_force_donor, name_file, output, "down", "force_donor")
    """
    # ------------------------ Creation of flaking intron size and force 5' figure for spliceosome exon list
    sf_list = list(group_factor.u1_factors) + ['DDX5_DDX17'] + ["CCE"]
    list_flanking_intron_size = []
    list_force_donor = []
    for sf_name in sf_list:
        if sf_name != exon_type:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx_sed, sf_name, "down")
            print("%s : %s exons" % (sf_name, len(exon_list)))
            new_exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list]
            list_force_donor.append(get_redundant_list_of_value(cnx_sed, new_exon_list, "force_donor"))
            upstream_intron_size = get_redundant_list_of_value(cnx_sed, new_exon_list, "upstream_intron_size", True)
            downstream_intron_size = get_redundant_list_of_value(cnx_sed, new_exon_list, "downstream_intron_size", True)
            flanking_intron = [np.nanmean(np.array([upstream_intron_size[i], downstream_intron_size[i]], dtype=float))
                              for i in range(len(upstream_intron_size))]
            flanking_intron = np.array(flanking_intron, dtype=float)
            flanking_intron = flanking_intron[~np.isnan(flanking_intron)]
            list_flanking_intron_size.append(flanking_intron)
        else:
            list_force_donor.append(get_control_exon_information(cnx_sed, exon_type, "force_donor"))
            upstream_intron_size = get_control_exon_information(cnx_sed, exon_type, "upstream_intron_size")
            downstream_intron_size = get_control_exon_information(cnx_sed, exon_type, "downstream_intron_size")
            print(upstream_intron_size[0:5])
            print(downstream_intron_size[0:5])
            flanking_intron = [np.nanmean(np.array([upstream_intron_size[i], downstream_intron_size[i]], dtype=float))
                               for i in range(len(upstream_intron_size))]
            flanking_intron = np.array(flanking_intron, dtype=float)
            flanking_intron = flanking_intron[~np.isnan(flanking_intron)]
            list_flanking_intron_size.append(flanking_intron)

    create_figure(list_force_donor, sf_list, output, "down", "force_donor_U1", "U1" )
    create_figure(list_flanking_intron_size, sf_list, output, "down", "flanking_intron_size_U1", "U1" )
    """
    # ------------------------ Creation of flaking intron size and force 3' and bp and ppt score figure for spliceosome U2 exon list
    sf_list = list(group_factor.u1_factors) + list(group_factor.u2_factors) + ["CCE"]
    list_flanking_intron_size = []
    list_force_acceptor = []
    dict_score_3ss = {i: {"bp_score_list": [], "ppt_score_list": [], "nb_bp_list": [], "nb_good_bp_list": []} for i in
                      [100, 50, 25]}
    for sf_name in sf_list:
        if sf_name != exon_type:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx_sed, sf_name, "down")
            print("%s : %s exons" % (sf_name, len(exon_list)))
            new_exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list]
            list_force_acceptor.append(get_redundant_list_of_value(cnx_sed, new_exon_list, "force_acceptor"))
            for size in dict_score_3ss.keys():
                bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list = get_bp_ppt_score_list(ctrl_output, new_exon_list, sf_name, size)
                dict_score_3ss[size]["bp_score_list"].append(bp_score_list)
                dict_score_3ss[size]["ppt_score_list"].append(ppt_score_list)
                dict_score_3ss[size]["nb_bp_list"].append(nb_bp_list)
                dict_score_3ss[size]["nb_good_bp_list"].append(nb_good_bp_list)
            upstream_intron_size = get_redundant_list_of_value(cnx_sed, new_exon_list, "upstream_intron_size", True)
            downstream_intron_size = get_redundant_list_of_value(cnx_sed, new_exon_list, "downstream_intron_size", True)
            flanking_intron = [np.nanmean(np.array([upstream_intron_size[i], downstream_intron_size[i]], dtype=float))
                              for i in range(len(upstream_intron_size))]
            flanking_intron = np.array(flanking_intron, dtype=float)
            flanking_intron = flanking_intron[~np.isnan(flanking_intron)]
            list_flanking_intron_size.append(flanking_intron)
        else:
            for size in dict_score_3ss.keys():
                mod = __import__("%s_%s_bp_ppt_score" % (exon_type, size))
                dict_score_3ss[size]["bp_score_list"].append(mod.bp_score)
                dict_score_3ss[size]["ppt_score_list"].append(mod.ppt_score)
                dict_score_3ss[size]["nb_bp_list"].append(mod.nb_bp)
                dict_score_3ss[size]["nb_good_bp_list"].append(mod.nb_good_bp)
            list_force_acceptor.append(get_control_exon_information(cnx_sed, exon_type, "force_acceptor"))
            upstream_intron_size = get_control_exon_information(cnx_sed, exon_type, "upstream_intron_size")
            downstream_intron_size = get_control_exon_information(cnx_sed, exon_type, "downstream_intron_size")
            flanking_intron = [np.nanmean(np.array([upstream_intron_size[i], downstream_intron_size[i]], dtype=float))
                               for i in range(len(upstream_intron_size))]
            flanking_intron = np.array(flanking_intron, dtype=float)
            flanking_intron = flanking_intron[~np.isnan(flanking_intron)]
            list_flanking_intron_size.append(flanking_intron)

    create_figure(list_force_acceptor, sf_list, ctrl_output, "down", "force_acceptor_U2", "U2" )
    create_figure(list_flanking_intron_size, sf_list, ctrl_output, "down", "flanking_intron_size_U2", "U2" )
    for size in dict_score_3ss.keys():
        # create_figure(dict_score_3ss[size]["ppt_score_list"], sf_list, output, "down", "ppt_score_U2(%snt)" % size, "U2")
        # create_figure(dict_score_3ss[size]["bp_score_list"], sf_list, output, "down", "bp_score_U2(%snt)" % size, "U2" )
        # create_figure(dict_score_3ss[size]["nb_bp_list"], sf_list, output, "down", "nb_branch_point_U2(%snt)" % size, "U2" )
        # create_figure(dict_score_3ss[size]["nb_good_bp_list"], sf_list, output, "down", "nb_good_branch_point_U2(%snt)" % size, "U2" )
        create_barplot(dict_score_3ss[size]["nb_bp_list"], sf_list, output, "down", "prop_nb_branch_point_spliceosome(%snt)" % size)
        create_barplot(dict_score_3ss[size]["nb_good_bp_list"], sf_list, output, "down", "prop_nb_good_branch_point_spliceosome(%snt)" % size)


if __name__ == "__main__":
    main()