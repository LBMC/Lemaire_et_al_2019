#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:

    The goal of this script is to generate correlation figure between the median frequency of each project \
    ant he gene_size or median intron size for each project only focusing on up, down and up and up and down \
    regulated exons.
"""

import figure_producer
import numpy as np
import os
import exon_control_handler
import plotly
import plotly.graph_objs as go
from scipy import stats
import argparse
import union_dataset_function
import math

def get_median_value(cnx, id_projects, target_column, control_dic, regulation, nt=None):
    """
    Return the median value of target_column in ``regulation`` exons of every  ``id_projects``.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param id_projects:  (list of int) the splicing lore id projects of every projects of interest
    :param target_column: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param control_dic: (dictionnary of list of float) median value of each possible control exons of \
    each feature in sed database.
    :param regulation: (list of string) up or down or up + down
    :param nt: (string) the nt of interest
    :return: (float) the relative median value (compared to control exons)  of ``target_column`` for every \
    ``regulation`` exons in every projects ``id_projects``
    """
    values_list = []
    for i in range(len(id_projects)):
        exon_list = []
        for j in range(len(regulation)):
            exon_list += figure_producer.get_ase_events(cnx, id_projects[i], regulation[j])
        if nt:
            values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt))
        else:
            values = np.array(figure_producer.get_list_of_value(cnx, exon_list, target_column))
        median_obs = np.median(values[~np.isnan(values)])
        if nt:
            final_value = float(median_obs - control_dic[target_column][nt]) / control_dic[target_column][nt] * 100
        else:
            #final_value = float(math.log(median_obs) - math.log(control_dic[target_column])) / math.log(control_dic[target_column]) * 100
            final_value = float(median_obs - control_dic[target_column]) / control_dic[target_column] * 100
        values_list.append(final_value)
    return values_list


def get_exons_values(cnx, sf_list, target_column, control_dic, regulation, nt):
    """
    Return the relative value (when compared to the median value of CCE exon) of target_column in `\
    `regulation`` exons regulated by a splicing factor in different cell lines.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param target_column: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param control_dic: (dictionnary of list of float) median value of each possible control exons of \
    each feature in sed database.
    :param regulation: (list of string) up or down or up + down
    :param nt: (string) the nt of interest
    :return: 3 lists :

        * values : (list of list of float) each sublist corresponds to the value of `` target_column`` for \
        every exons regulated by a splicing factor
        * exon_name : (list of list of string) each sublist corresponds to the name of \
        every exons regulated by a splicing factor - the value in the sublist **i** position **j** \
         in the ``value`` and ``exon_name`` corresponds to the same exons
        * all_sf (list of string) list of each sf studied
    """
    all_sf = []
    exon_name = []
    values = []
    for sf_name in sf_list:
        exon_list = []
        for j in range(len(regulation)):
            exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation[j])
        exon_name.append(["%s_%s" % (union_dataset_function.get_gene_name(cnx, a[0]), a[1]) for a in exon_list])
        all_sf.append(sf_name)
        if nt:
            values.append(np.array(figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt)))
        else:
            values.append(np.array(figure_producer.get_redundant_list_of_value(cnx, exon_list, target_column)))
    return values, exon_name, all_sf


def color_maker(name_projects, value_size, value_iupac):
    """
    Split the data in 4 group according to the splicing factor they regulate.

    :param name_projects: (list of string) the name of the project of interest
    :param value_size: (list of float) the list of relative median value of ``size_scale`` for \
    every projects named ``projects_names``.
    :param value_iupac: (list of float) the list of relative median value of the frequency in ``nt`` for \
    every projects named ``projects_names``.
    :return: (dictionary of 4 list of float or string) for each group (key of the dictionary) we have, \
     the x coordinates of the exons group, the y corrdinates, the color code of the group, and the \
     name of every project in the group.
    """
    spliceosome = ["PRPF6", "PRPF8", "SF1", "SF3A3", "SF3B1", "SF3B4", "SNRNP200", "SNRNP40",
                   "SRNRP70", "SNRPC", "U2AF1", "U2AF2"]
    result = {'HNRNP': [[], [], ['rgba(45, 165, 43, 0.8)'], []],
              'Misc': [[], [], ['rgba(185, 73, 184, 0.8)'], []],
              'SRSF': [[], [], ['rgba(215, 126, 0, 0.8)'], []],
              'RBM': [[], [], ['rgba(247, 80, 70, 0.8)'], []],
              'Spliceosome': [[], [], ['rgba(55, 126, 219, 0.8)'], []]}
    for i in range(len(name_projects)):
        if "srsf" in name_projects[i].lower():
            result["SRSF"][0].append(value_size[i])
            result["SRSF"][1].append(value_iupac[i])
            result["SRSF"][3].append(name_projects[i])
        elif "hnrn" in name_projects[i].lower():
            result["HNRNP"][0].append(value_size[i])
            result["HNRNP"][1].append(value_iupac[i])
            result["HNRNP"][3].append(name_projects[i])
        elif "rbm" in name_projects[i].lower():
            result["RBM"][0].append(value_size[i])
            result["RBM"][1].append(value_iupac[i])
            result["RBM"][3].append(name_projects[i])
        else:
            test = "ko"
            for sp_name in spliceosome:
                if sp_name in name_projects[i].upper() and test == "ko":
                    result["Spliceosome"][0].append(value_size[i])
                    result["Spliceosome"][1].append(value_iupac[i])
                    result["Spliceosome"][3].append(name_projects[i])
                    test = "ok"
            if test == "ko":  # no spliceosome componant
                result["Misc"][0].append(value_size[i])
                result["Misc"][1].append(value_iupac[i])
                result["Misc"][3].append(name_projects[i])
    return result


def figure_creator(values_size, values_iupac, projects_names, regulation, size_scale, nt_name, ctrl, output, type=None):
    """
    Create a scatter plot showing the potential correlation between projects.

    :param values_size: (list of float) the list of relative median value of ``size_scale`` for \
    every projects named ``projects_names``.
    :param values_iupac: (list of float) the list of relative median value of the frequency in ``nt`` for \
    every projects named ``projects_names``.
    :param projects_names: (list of string) list of every projects name studied.
    :param regulation: (list of string) the regulation chosen up down or up and down
    :param size_scale: (string) either gene_size or median_intron_size.
    :param nt_name: (string) the name of the nucleotide of interest
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    :param type: (string or NoneType object) the type of legend to put in the graphics
    """
    trace_pattern = color_maker(projects_names, values_size, values_iupac)
    data = []
    slope, intercept, r_value, p_value, std_err = stats.linregress(values_size, values_iupac)
    line = slope * np.array(values_size) + intercept
    p2 = go.Scatter(x=values_size,
                    y=line,
                    mode='lines',
                    line=dict(color='red', width=3),
                    name="Fit"
                    )
    data.append(p2)
    for key in trace_pattern:
        data.append(go.Scatter(
            x=trace_pattern[key][0],
            y=trace_pattern[key][1],
            name=key,
            mode='markers',
            text=trace_pattern[key][3],
            marker=dict(
                size=10,
                color=trace_pattern[key][2][0],
                line=dict(width=1))
        ))
    cor, pval = stats.pearsonr(values_size, values_iupac)

    if not type:
        main_title = 'Correlation between %s and %s frequency in genes containing %s-regulated exons in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (size_scale, nt_name, regulation, ctrl, round(cor, 2), pval)
        x_title = 'relative %s' % size_scale
        y_title = 'relative %s frequency' % nt_name
        figname = '%s%s_%s_correlation_graphs_%s.html'% (output, nt_name, size_scale, regulation[0])
    elif type == 1:
        main_title = 'Correlation between down %s frequency and up %s frequency for exons in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (size_scale, nt_name, ctrl, round(cor, 2), pval)
        x_title = 'relative %s frequency - down exons' % size_scale
        y_title = 'relative %s frequency - up exons' % nt_name
        figname = '%s%s_correlation_graphs.html'% (output, nt_name)
    elif type == 2:
        main_title = 'Correlation between gene %s nt frequency and exon %s nt frequency for %s exons in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (size_scale, nt_name, regulation[0], ctrl, round(cor, 2), pval)
        x_title = 'relative %s frequency - %s gene' %(size_scale, regulation[0])
        y_title = 'relative %s frequency - %s exons' % (nt_name, regulation[0])
        figname = '%s%s_gene_vs_exon_correlation_graphs_%s.html'% (output, nt_name, regulation[0])
    elif type == 3:
        main_title = 'Correlation between exon %s nt frequency and %s %s nt frequency for %s exons in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (
        nt_name, size_scale, nt_name, regulation[0], ctrl, round(cor, 2), pval)
        x_title = 'relative %s frequency in exon - %s exons' % (nt_name, regulation[0])
        y_title = 'relative %s frequency in %s - %s exons' % (nt_name, size_scale, regulation[0])
        figname = '%s%s_intron_vs_exon_correlation_graphs_%s.html' % (output, nt_name, regulation[0])
    else:
        main_title = 'Correlation between the relative %s and the %s nt frequency in %s exons for every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (
        size_scale, nt_name, regulation[0], ctrl, round(cor, 2), pval)
        x_title = 'relative %s - %s exons' % (size_scale, regulation[0])
        y_title = 'relative %s frequency in exons - %s exons' % (nt_name, regulation[0])
        figname = '%s%s_%s_vs_exon_correlation_graphs_%s.html' % (output, nt_name, size_scale, regulation[0])
    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title),
        yaxis=dict(
            title=y_title),
        showlegend=True
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname,
                        auto_open=False)


def figure_creator_exon(values_feature1, values_iupac, all_sf, regulation, size_scale, nt_name, ctrl, exon_name, output, type=None):
    """
    Create a scatter plot showing the potential correlation for every exons regulated each splicing factor in \
    multiple cell lines.

    :param values_feature1: (list of float) the list of relative value of a particular feature for \
    exons ``regulation`` regulated by ``sf_name``
    :param values_iupac: (list of float) the list of relative value of the frequency in ``nt`` for \
    exons ``regulation`` regulated by ``sf_name``
    :param sf_name: (list of string) list of every splicing factor studied.
    :param regulation: (list of string) the regulation chosen up down or up and down
    :param size_scale: (string) either gene_size or median_intron_size.
    :param nt_name: (string) the name of the nucleotide of interest
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    :param exon_name: (list of string) the name of the interest exons
    :param type: (string or NoneType object) the type of legend to put in the graphics
    """
    log=""
    if max(values_feature1[0]) >= 1000:
        values_feature1 = [list(map(math.log10, listval)) for listval in values_feature1]
        log = "log10"
    data = []
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(all_sf))]
    for i in range(len(all_sf)):
        slope, intercept, r_value, p_value, std_err = stats.linregress(values_feature1[i], values_iupac[i])
        line = slope * np.array(values_feature1[i]) + intercept
        cor, pval = stats.pearsonr(values_feature1[i], values_iupac[i])
        data.append(go.Scatter(x=values_feature1[i],
                        y=line,
                        visible="legendonly",
                        mode='lines',
                        line=dict(color=c[i], width=3),
                        name="Fit %s p=%.2E, c=%.2f" % (all_sf[i], pval, cor)
                        ))
        data.append(go.Scatter(
            x=values_feature1[i],
            y=values_iupac[i],
            name="%s (%s exons)" % (all_sf[i], len(values_feature1[i])),
            mode='markers',
            visible="legendonly",
            text=exon_name[i],
            marker=dict(
                size=10,
                color=c[i],
                line=dict(width=1))
        ))


    if not type:
        main_title = 'Correlation between %s and %s frequency in genes containing %s-regulated (relative value against %s control - EXON LEVEL)' % (size_scale, nt_name, regulation, ctrl)
        x_title = '%s %s' % (log, size_scale)
        y_title = '%s frequency' % nt_name
        figname = '%s%s_%s_correlation_graphs_%s.html'% (output, nt_name, size_scale, regulation)
    elif type == 2:
        main_title = 'Correlation between gene %s nt frequency and exon %s nt frequency for %s exons in every splicing lore project<br> (relative value against %s control - EXON LEVEL)' % (size_scale, nt_name, regulation[0], ctrl)
        x_title = '%s %s frequency - %s gene' % (log, size_scale, regulation[0])
        y_title = '%s frequency - %s exons' % (nt_name, regulation[0])
        figname = '%s%s_gene_vs_exon_correlation_graphs_%s.html'% (output, nt_name, regulation[0])
    elif type == 3:
        main_title = 'Correlation between exon %s nt frequency and %s %s nt frequency for %s exons in every splicing lore project<br> (relative value against %s control - EXON LEVEL)' % (
        nt_name, size_scale, nt_name, regulation[0], ctrl)
        x_title = '%s %s frequency in exon - %s exons' % (log, nt_name, regulation[0])
        y_title = '%s frequency in %s - %s exons' % (nt_name, size_scale, regulation[0])
        figname = '%s%s_intron_vs_exon_correlation_graphs_%s.html' % (output, nt_name, regulation[0])
    else:
        main_title = 'Correlation between the relative %s and the %s nt frequency in %s exons for every splicing lore project<br> (relative value against %s control - EXON LEVEL)' % (
        size_scale, nt_name, regulation[0], ctrl)
        x_title = '%s %s - %s exons' % (log, size_scale, regulation[0])
        y_title = '%s frequency in exons - %s exons' % (nt_name, regulation[0])
        figname = '%s%s_%s_vs_exon_correlation_graphs_%s.html' % (output, nt_name, size_scale, regulation[0])
    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title),
        yaxis=dict(
            title=y_title),
        showlegend=True
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname,
                        auto_open=False)


def remove_none_values(values_list1, values_list2, list_exon_name):
    """
    Remove the none values of each sublist of this parameters.

    :param values_list1: (list of list of float) each float value of a sublist corresponds to a characteristic \
    of a particular exon regulated by a splicing factor
    :param values_list2:  (list of list of float) each float value of a sublist corresponds to a characteristic \
    of a particular exon regulated by a splicing factor
    :param list_exon_name: (list of list of string) each value of a substring si an exon name
    :return: the list without the none values
    """
    values_list1 = [list(a) for a in values_list1]
    values_list2 = [list(a) for a in values_list2]
    list_exon_name = [list(a) for a in list_exon_name]
    i = 0
    while i < len(values_list1):
        j = 0
        while j < len(values_list1[i]):
            if values_list1[i][j] is None or values_list2[i][j] is None:
                del(values_list1[i][j])
                del(values_list2[i][j])
                del(list_exon_name[i][j])
                j -= 1
            j += 1
        i += 1
    return values_list1, values_list2, list_exon_name


def main(exon_level):
    """
    Create the correlation matrix (gene_size vs iupac)

    :param exon_level: (boolean) True to work on exon level False else
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    id_projects, name_projects = figure_producer.get_interest_project(cnx)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    regulations = [["up"], ["down"], ["up", "down"]]
    targets1 = ["gene_size", "median_intron_size"]
    target2 = "iupac_gene"
    if not exon_level:
        output_niv0 = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_size_frequency/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output_niv0):
            os.mkdir(output_niv0)
        for my_regulation in regulations:
            output_niv1 = output_niv0 + my_regulation[0] + "/"
            if len(my_regulation) > 1:
                name_reg = ["-".join(my_regulation)]
                output_niv1 = output_niv0 + name_reg[0] + "/"
            if not os.path.isdir(output_niv1):
                os.mkdir(output_niv1)
            for target1 in targets1:
                value_target1 = get_median_value(cnx, id_projects, target1, ctrl_dic, my_regulation, nt=None)
                for nt in nt_list:
                    value_target2 = get_median_value(cnx, id_projects, target2, ctrl_dic, my_regulation, nt=nt)
                    if len(my_regulation) > 1:
                        figure_creator(value_target1, value_target2, name_projects, name_reg[0],
                                       target1, nt, exon_type, output_niv1)
                    else:
                        figure_creator(value_target1, value_target2, name_projects, my_regulation[0],
                                       target1, nt, exon_type, output_niv1)
    else:
        output_niv0 = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_size_frequency_exon_level/"
        list_sf = sorted(union_dataset_function.get_splicing_factor_name(cnx))
        if not os.path.isdir(output_niv0):
            os.mkdir(output_niv0)
        for my_regulation in regulations:
            output_niv1 = output_niv0 + my_regulation[0] + "/"
            if len(my_regulation) > 1:
                name_reg = ["-".join(my_regulation)]
                output_niv1 = output_niv0 + name_reg[0] + "/"
            if not os.path.isdir(output_niv1):
                os.mkdir(output_niv1)
            for target1 in targets1:
                value_target1, exon_name, all_sf = get_exons_values(cnx, list_sf, target1, ctrl_dic, my_regulation, nt=None)
                for nt in nt_list:
                    value_target2, exon_name, all_sf = get_exons_values(cnx, list_sf, target2, ctrl_dic, my_regulation, nt=nt)
                    if len(my_regulation) > 1:
                        figure_creator_exon(value_target1, value_target2, all_sf, name_reg[0],
                                       target1, nt, exon_type, exon_name, output_niv1)
                    else:
                        figure_creator_exon(value_target1, value_target2, all_sf, my_regulation[0],
                                       target1, nt, exon_type, exon_name, output_niv1)


def main_up_vs_down():
    """
    Create the correlation matrix (iupac up vs iupac down).
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    id_projects, name_projects = figure_producer.get_interest_project(cnx)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_up_vs_down/"
    # If the output directory does not exist, then we create it !
    if not os.path.isdir(output):
        os.mkdir(output)
    for nt in nt_list:
        value_target2 = get_median_value(cnx, id_projects, "iupac_exon", ctrl_dic, ["up"], nt=nt)
        value_target1 = get_median_value(cnx, id_projects, "iupac_exon", ctrl_dic, ["down"], nt=nt)
        figure_creator(value_target1, value_target2, name_projects, None,
                       nt , nt, exon_type, output, 1)


def iupac_gene_vs_iupac_exon(exon_level):
    """
    Correlation graphics between the relative median iupac gene frequency and the relative median iupac exon frequency \
    for up or down-regulated exon in every splicing lore project.

    :param exon_level: (boolean) true if we want to work at the exon level, fasle else.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    regulations = [["up"], ["down"]]
    if not exon_level:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_iupac_gene_exon/"
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for my_regulation in regulations:
            for nt in nt_list:
                value_target1 = get_median_value(cnx, id_projects, "iupac_gene", ctrl_dic, my_regulation, nt=nt)
                value_target2 = get_median_value(cnx, id_projects, "iupac_exon", ctrl_dic, my_regulation, nt=nt)
                figure_creator(value_target1, value_target2, name_projects, my_regulation,
                               nt, nt, exon_type, output, 2)
    else:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_iupac_gene_exon_EXON_LEVEL/"
        list_sf = sorted(union_dataset_function.get_splicing_factor_name(cnx))
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for my_regulation in regulations:
            for nt in nt_list:
                value_target1, exon_name, all_sf = get_exons_values(cnx, list_sf, "iupac_gene",
                                                                    ctrl_dic, my_regulation, nt=nt)
                value_target2, exon_name, all_sf = get_exons_values(cnx, list_sf, "iupac_exon",
                                                                    ctrl_dic, my_regulation, nt=nt)
                figure_creator_exon(value_target1, value_target2, all_sf, my_regulation,
                                    nt, nt, exon_type, exon_name, output, 2)

def iupac_exon_vs_iupac_intron_proxi(exon_level):
    """
    Correlation graphics between the relative median iupac exon frequency and the relative median iupac intron \
    upstream and downstream frequency \
    for up or down-regulated exon in every splicing lore project.s

    :param exon_level: (boolean) true for exon level (see every exon in graphics) false else.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    regulations = [["up"], ["down"]]
    intron_targets = ["iupac_upstream_intron_proxi", "iupac_downstream_intron_proxi"]
    if not exon_level:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_iupac_exon_intron_proxi/"
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for my_intron_target in intron_targets:
            output_final = output +  my_intron_target + "/"
            if not os.path.isdir(output_final):
                os.mkdir(output_final)
            for my_regulation in regulations:
                for nt in nt_list:
                    value_target1 = get_median_value(cnx, id_projects, "iupac_exon", ctrl_dic, my_regulation, nt=nt)
                    value_target2 = get_median_value(cnx, id_projects, my_intron_target, ctrl_dic, my_regulation, nt=nt)
                    figure_creator(value_target1, value_target2, name_projects, my_regulation,
                                   my_intron_target , nt, exon_type, output_final, 3)
    else:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_iupac_exon_intron_proxi_EXON_LEVEL/"
        list_sf = sorted(union_dataset_function.get_splicing_factor_name(cnx))
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for my_intron_target in intron_targets:
            output_final = output + my_intron_target + "/"
            if not os.path.isdir(output_final):
                os.mkdir(output_final)
            for my_regulation in regulations:
                for nt in nt_list:
                    value_target1, exon_name, all_sf = get_exons_values(cnx, list_sf, "iupac_exon",
                                                                        ctrl_dic, my_regulation, nt=nt)
                    value_target2, exon_name, all_sf = get_exons_values(cnx, list_sf, my_intron_target,
                                                                        ctrl_dic, my_regulation, nt=nt)
                    value_target1, value_target2, exon_name = remove_none_values(value_target1, value_target2, exon_name)
                    figure_creator_exon(value_target1, value_target2, all_sf, my_regulation,
                                        my_intron_target, nt, exon_type, exon_name, output_final, 3)

def  force_vs_iupac_exon(exon_level):
    """
    Correlation graphics between the 5' and 3' force and the relative median iupac exon frequency \
    for up or down-regulated exon in every splicing lore project.

    :param exon_level: (boolean) true if you want to work at the exon level, false else.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    regulations = [["up"], ["down"]]
    forces = ["force_acceptor", "force_donor"]
    if not exon_level:
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_force_vs_iupac_exon/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for i in range(len(forces)):
            my_force = forces[i]
            output_final = output +  forces[i] + "/"
            if not os.path.isdir(output_final):
                os.mkdir(output_final)
            for my_regulation in regulations:
                for nt in nt_list:
                    value_target1 = get_median_value(cnx, id_projects, my_force, ctrl_dic, my_regulation, nt=None)
                    value_target2 = get_median_value(cnx, id_projects, "iupac_exon", ctrl_dic, my_regulation, nt=nt)
                    figure_creator(value_target1, value_target2, name_projects, my_regulation,
                                   my_force , nt, exon_type, output_final, 4)
    else:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_force_vs_iupac_exon_EXON_LEVEL/"
        list_sf = sorted(union_dataset_function.get_splicing_factor_name(cnx))
        if not os.path.isdir(output):
            os.mkdir(output)
        for i in range(len(forces)):
            my_force = forces[i]
            output_final = output +  forces[i] + "/"
            if not os.path.isdir(output_final):
                os.mkdir(output_final)
            for my_regulation in regulations:
                for nt in nt_list:
                    value_target1, exon_name, all_sf = get_exons_values(cnx, list_sf, my_force,
                                                                        ctrl_dic, my_regulation, nt=None)
                    value_target2, exon_name, all_sf = get_exons_values(cnx, list_sf, "iupac_exon",
                                                                        ctrl_dic, my_regulation, nt=nt)
                    value_target1, value_target2, exon_name = remove_none_values(value_target1, value_target2, exon_name)
                    figure_creator_exon(value_target1, value_target2, all_sf, my_regulation,
                                        my_force , nt, exon_type, exon_name, output_final, 4)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
   The goal of this script is to create correlation figure.
   The figures can show the correlation between:
        1 - The median_intron_size (and gene_size) and the iupac frequencies for every splicing lore project \
        (ech dot corresponds to the median iupac frequency and the median intron size of every up or down exons in \
        a splicing lore project)
        2 - The iupac frequency between up and down regulated exon for every splicing lore project. (A dot corresponds \
        to  the median iupac frequency for up and down exons

    """)
    # Arguments for the parser

    parser.add_argument('--exon_level', dest="exon_level", help="True if you want to create correlation at the"
                                                                "exon  level, false else", default=False)
    req_arg = parser.add_argument_group("required arguments")

    req_arg.add_argument('--type', '-t', dest='type', help="the type of graphic you want to create (gene_level' or"
                                                           "up_vs_down",
                         required=True)
    args = parser.parse_args()

    if args.exon_level == "False":
        args.exon_level = False
    if args.exon_level == "True":
        args.exon_level = True


    my_types = ["gene_size_vs_gene_iupac", "iupac_up_vs_iupac_down", "iupac_gene_vs_iupac_exon", "iupac_exon_vs_iupac_intron_proxi", "force_vs_iupac_exon"]
    if args.type not in my_types:
        parser.error("Wrong parameter type\n It can only be %s" % " or ".join(my_types))
    if args.type == "gene_size_vs_gene_iupac":
        main(args.exon_level)
    elif args.type == "iupac_up_vs_iupac_down":
        main_up_vs_down()
    elif args.type == "iupac_gene_vs_iupac_exon":
        iupac_gene_vs_iupac_exon(args.exon_level)
    elif args.type == "iupac_exon_vs_iupac_intron_proxi":
        iupac_exon_vs_iupac_intron_proxi(args.exon_level)
    elif args.type == "force_vs_iupac_exon":
        force_vs_iupac_exon(args.exon_level)

if __name__ == "__main__":
    launcher()