#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:

    This program will see if a list of exons GC rich is more often regulated by u1 than AT-rich or CCE:ACE... list\
    of exons.
"""

import sys
sys.path.insert(0, "/media/nicolas/DD_1/Splicing_Lore_project/GC_rich_AT_rich_exon_list/src")
import group_factor
import numpy as np
import os
import sqlite3
import union_dataset_function
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import plotly.figure_factory as ff
import plotly
import plotly.graph_objs as go

def get_exons_list(cnx, sf_list, regulation):
    """
    Return every non-redundant exons regulated by  at least one factor in ``sf_list`` (with the regulation \
    ``regulation``)

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param regulation: (list of string) up or down or up + down
    :return: (list of list of int) list of exons shownig the regulation ``regulation`` at least for a factor \
    in ``sf_list``
    """
    exon_list = []
    for sf_name in sf_list:
        exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_list = union_dataset_function.washing_events_all(exon_list)
    return ["_".join(map(str, b)) for b in exon_list]


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


def exon_difference(exon_list1, exon_list2):
    """
    Return the difference of 2 exon lists

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the difference between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon not in exon_list2]
    return exon_list


def get_control_exon(cnx, exon_type):
    """
    Get the wanted control exons
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return:
        * result: (list of 2 int) list of control exons
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos
                   FROM sed"""
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    values = []
    for val in result:
        values.append("_".join(list(map(str, list(val)))))
    return values


def simplify_dic(dic):
    """
    For each key of ``dic`` get the mean value of the list of values linked to this key.
    :param dic: (dic of list of float) dic of deltapsi (values) of exons (key)
    :return: (dic of float) dic of average deltapsi of exons.
    """
    new_dic = {}
    for key in dic.keys():
       new_dic[key] = np.mean(dic[key])
    return new_dic


def get_dpsi_of_list_of_factor(cnx, sf_list):
    """
    get the delta_psi of a list of factors
    :param cnx: (sqlite3 connect object) connection to seddb
    :param sf_list: (list of string) list of splicing factors
    :return: (dictionary of float) link for each exon their delta_psi
    """
    dpsi_dic_list = {}
    list_id = []
    cursor = cnx.cursor()
    for sf_name in sf_list:
        list_id += union_dataset_function.get_projects_links_to_a_splicing_factor(cnx, sf_name)
    for id_project in  list_id:
        query = """SELECT gene_id, exon_skipped, delta_psi
                   FROM ase_event
                   WHERE id_project = {}
                   AND delta_psi IS NOT NULL""".format(id_project)
        cursor.execute(query)
        result = cursor.fetchall()
        print(id_project, ":", len(result))
        for exon in result:
            exon_name = "_".join(list(map(str, exon[0:2])))
            if exon_name not in dpsi_dic_list.keys():
                dpsi_dic_list[exon_name] = [float(exon[2])]
            else:
                dpsi_dic_list[exon_name].append(float(exon[2]))
    return simplify_dic(dpsi_dic_list)


def mann_withney_test_r(gc_content_test, gc_control):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(gc_content_test), v.FloatVector(gc_control))[0])
    return pval


def make_dpsi_density_hist(dpsi_list_1, dpsi_list_2, name_list1, name_list2, name_ctrl, output):
    """
    Create a density plot showing the deltapsi in ``name_ctrl`` of exon exon _list ``name_list1`` \
    and ``name_list2``
    :param dpsi_list_1: (list of float) list of delta_psi in ``name_ctrl`` of exon ``name_list1``
    :param dpsi_list_2:(list of float) list of delta_psi in ``name_ctrl`` of exon ``name_list2``
    :param name_list1: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_1`` \
    in ``name_ctrl``
    :param name_list2: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_2`` \
    in ``name_ctrl``
    :param name_ctrl: (string) the project for which we have the dpsi ``dpsi_list_1`` and ``dpsi_list_2``
    :param output: (string) path where figure will be created
    """
    pval = mann_withney_test_r(dpsi_list_1, dpsi_list_2)
    hist_data = [dpsi_list_1, dpsi_list_2]
    group_labels = [name_list1, name_list2]
    fig = ff.create_distplot(hist_data, group_labels,
                             bin_size=2.5, curve_type='kde', show_rug =False, show_hist=False,
                             colors=['#1f77b4', '#2ca02c'])
    main_title = 'Distribution of deltapsi in %s between %s and %s exon_list<br> wilcoxon rank test pvalue = %s' % (name_ctrl, name_list1, name_list2, pval)

    x_title = 'deltapsi'
    y_title = 'Frequency'
    figname = '%sdensity_%s-dpsi_%s_exon_vs_%s.html' % (output, name_ctrl, name_list1, name_list2)
    fig['layout'].update(title=main_title)
    fig['layout'].update(xaxis=dict(title=x_title))
    fig['layout'].update(yaxis=dict(title=y_title))
    plotly.offline.plot(fig, filename=figname, auto_open=False)


def make_dpsi_density_hist_ctrl(val_mean_dpsi, ctrl_pval, name_list1, name_list2, name_ctrl, output, nb_iteration):
    """
    Create a density plot showing the deltapsi in ``name_ctrl`` of exon exon _list ``name_list1`` \
    and ``name_list2``
    :param val_mean_dpsi: (float) the mean value of dpsi for exon regulatedby ``name_list1`` (dpsi value taken \
    from ``name_ctrl``
    :param ctrl_pval:(list of float) list of delta_psi in ``name_ctrl`` of ctrl exon ``name_list2``
    in ``name_ctrl``
    :param name_list1: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_1`` \
    in ``name_ctrl``
    :param name_list2: (string) the name of the list of exon that have the list of delta_psi ``dpsi_list_2`` \
    in ``name_ctrl``
    :param name_ctrl: (string) the project for which we have the dpsi ``dpsi_list_1`` and ``dpsi_list_2``
    :param output: (string) path where figure will be created
    :param nb_iteration: (int) the number of iteration
    """
    pval = min(sum(np.array(ctrl_pval) >= val_mean_dpsi) / len(ctrl_pval),
               sum(np.array(ctrl_pval) <= val_mean_dpsi) / len(ctrl_pval))
    hist_data = [ctrl_pval]
    group_labels = [name_list2]
    fig = ff.create_distplot(hist_data, group_labels,
                             curve_type='kde', show_rug =False, show_hist=False,
                             colors=['#1f77b4', '#2ca02c'])
    main_title = 'Distribution of deltapsi in %s between %s and %s exon_list<br> iteration %s - pvalue = %s' % (name_ctrl, name_list1, name_list2, nb_iteration, pval)

    x_title = 'deltapsi'
    y_title = 'Frequency'
    figname = '%sdensity_%s-dpsi_%s_exon_vs_%s.html' % (output, name_ctrl, name_list1, name_list2)
    fig['layout'].update(title=main_title)
    fig['layout'].update(xaxis=dict(title=x_title))
    fig['layout'].update(yaxis=dict(title=y_title))
    fig['layout'].update(shapes=[dict(type="line", x0=val_mean_dpsi, y0=0, x1=val_mean_dpsi, y1=30, line=dict(color="red", width=2))])
    plotly.offline.plot(fig, filename=figname, auto_open=False)


def analysis_maker(gc_exon_list, at_exon_list, gc_pure_exon_list, at_pure_exon_list, splicesome_dic,
                   nb_iteration, output):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list, one containing exons regulated by splicing factors regulating AT rich exons and the other \
    containing exons regulated by splicing factors regulating GC rich exons
    :param gc_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating GC rich exons.
    :param at_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating AT rich exons.
    :param gc_pure_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating \
    GC rich exons and not regulated by splicing factors regulating AT rich exons.
    :param at_pure_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating \
    AT rich exons and not regulated by splicing factors regulating GC rich exons.
    :param splicesome_dic: (dict of dic of float) dictionary that contains every dic of delta psi exons regulated by a \
    particular U1 spliceosome factors.
    :param nb_iteration: (int) the number of iteration we gonna make
    :param output: (string) path where the result fig will be created
    """
    output += "GC_vs_AT_fig/"
    if not os.path.isdir(output):
        os.mkdir(output)
    at_gc = ["pure", "all"]
    ctrl = ["all"]
    for gc_at_set in at_gc:
        if gc_at_set == "all":
            cur_gc_list = gc_exon_list
            cur_at_list = at_exon_list
        else:
            cur_gc_list = gc_pure_exon_list
            cur_at_list = at_pure_exon_list
        for ctrl_set in ctrl:
            for key in splicesome_dic.keys():
                if ctrl_set == "all":
                    cur_ctrl_keys = splicesome_dic[key].keys()
                list_gc_dpsi = [splicesome_dic[key][exon] for exon in cur_gc_list if exon in cur_ctrl_keys]
                list_at_dpsi = [splicesome_dic[key][exon] for exon in cur_at_list if exon in cur_ctrl_keys]
                make_dpsi_density_hist(list_gc_dpsi, list_at_dpsi, "gc_%s_group" % gc_at_set,
                                       "at_%s_group" % gc_at_set, "%s_%s" % (key, ctrl_set), output)


def analysis_maker_ctrl(gc_exon_list, gc_pure_exon_list, exon_ctrl, test_set, exon_type, splicesome_dic,
                   nb_iteration, output):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list, one containing exons regulated by splicing factors regulating AT rich exons and the other \
    containing exons regulated by splicing factors regulating GC rich exons
    :param gc_exon_list: (list of string) list of exons regulated by splicing factors regulating GC rich exons..
    :param gc_pure_exon_list: (list of string) list of exons regulated by splicing factors regulating \
    GC rich exons and not regulated by splicing factors regulating AT rich exons.
    :param exon_ctrl: (list of string) list of exons
    :param test_set: (string) the type of test_set used
    :param exon_type: (string) list of exons of interest
    :param splicesome_dic: (dict of dic of float) dictionary that contains every dic of delta psi exons regulated by a \
    particular U1 spliceosome factors.
    :param nb_iteration: (int) the number of iteration we gonna make
    :param output: (string) path where the result fig will be created
    """
    output += "%s_vs_%s_fig/" % (test_set, exon_type)
    if not os.path.isdir(output):
        os.mkdir(output)
    at_gc = ["pure", "all"]
    ctrl = ["all"]
    for gc_at_set in at_gc:
        if gc_at_set == "all":
            cur_gc_list = gc_exon_list
        else:
            cur_gc_list = gc_pure_exon_list
        for ctrl_set in ctrl:
            for key in splicesome_dic.keys():
                if ctrl_set == "all":
                    cur_ctrl_keys = splicesome_dic[key].keys()
                values_dpsi = list(splicesome_dic[key].values())
                mean_gc_dpsi = np.mean([splicesome_dic[key][exon] for exon in cur_gc_list if exon in cur_ctrl_keys])
                print("------", len(cur_gc_list), "--", len([splicesome_dic[key][exon] for exon in cur_gc_list if exon in cur_ctrl_keys]))
                mean_ctrl_dpsi = []
                print("-------------key : %s %s %s" % (key, gc_at_set, ctrl_set))
                for i in range(nb_iteration):
                    sys.stdout.write("iteration %s/%s     \r" % (i, nb_iteration))
                    list_ctrl_dpsi = np.random.choice(values_dpsi, len(cur_gc_list), replace=False)
                    mean_ctrl_dpsi.append(np.mean(list_ctrl_dpsi))
                # print(mean_ctrl_dpsi)
                make_dpsi_density_hist_ctrl(mean_gc_dpsi, mean_ctrl_dpsi, "%s_%s_group" % (test_set, gc_at_set),
                                       "%s_group" % exon_type, "%s_%s" % (key, ctrl_set), output, nb_iteration)


def main(test_set, comparison):
    """
    Make the enrichment analysis comparing the frequencies of exon regulated by splicesome factors \
     for an AT and GC exons list.

    :param comparison: (string) the comparison set of exon used
    :param test_set: (string) the test set used
    """
    nb_iteration = 10000
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2",
                                                                 "result/GC_AT_group_regulated_U1_U2/")
    if not os.path.isdir(output):
        os.mkdir(output)
    if test_set == "GC":
        div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
                     "U1" : group_factor.u1_factors, "U1_DDX5": list(group_factor.u1_factors) + ["DDX5_DDX17"],
                     "SNRPC": ["SNRPC"], "SNRNP70": ["SNRNP70"], "DDX5_17": ["DDX5_DDX17"]}
    else:
        div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
                     "U2" : group_factor.u2_factors, "SF1": ["SF1"], "U2AF1": ["U2AF1"], "U2AF2": ['U2AF2'],
                     "SF3B1": ["SF3B1"], "SF3B4": ["SF3B4"]}
    dic_exon = {}
    for name_group in ["AT_rich", "GC_rich"]:
        print("Getting all exon regulated by %s factor" % name_group)
        dic_exon[name_group] = get_exons_list(cnx, div_group[name_group], "down")
    at_gc_intersection = exon_intersection(dic_exon["AT_rich"], dic_exon["GC_rich"])
    # u1_u2_intersection = exon_intersection(get_exons_list(cnx, group_factor.u1_factors, "down"), get_exons_list(cnx, group_factor.u2_factors, "down"))
    dic_exon["GC_pure"] = exon_difference(dic_exon["GC_rich"], at_gc_intersection)
    dic_exon["AT_pure"] = exon_difference(dic_exon["AT_rich"], at_gc_intersection)
    for key in dic_exon:
        print("%s : %s" % (key, len(dic_exon[key])))
    dic_spliceosome = {}
    for key in div_group.keys():
        if "AT" not in key and "GC" not in key:
            dic_spliceosome[key] = get_dpsi_of_list_of_factor(cnx, div_group[key])
            print(key, ":", len(dic_spliceosome[key].keys()), "dpsi")
    if comparison == "AT" or comparison == "GC":
        analysis_maker(dic_exon["GC_rich"], dic_exon["AT_rich"], dic_exon["GC_pure"], dic_exon["AT_pure"], dic_spliceosome, nb_iteration, output)
    else:
        ctrl_exons = get_control_exon(cnx, comparison)
        if test_set == "GC":
            analysis_maker_ctrl(dic_exon["GC_rich"], dic_exon["GC_pure"], ctrl_exons, test_set, comparison, dic_spliceosome,
                   nb_iteration, output)
        else:
            analysis_maker_ctrl(dic_exon["AT_rich"], dic_exon["AT_pure"], ctrl_exons, test_set, comparison, dic_spliceosome,
                                nb_iteration, output)
    cnx.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

