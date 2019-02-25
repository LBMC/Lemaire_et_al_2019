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
import control_exon_adapter
import plotly
import plotly.graph_objs as go
from scipy import stats
import argparse
import union_dataset_function
import group_factor
import math
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7, "K": 8, "M": 9}
dnt_dic = {"AA": 0, "AC": 1, "AG": 2, "AT": 3, "CA": 4, "CC": 5,
           "CG": 6, "CT": 7, "GA": 8, "GC": 9, "GG": 10, "GT": 11,
           "TA": 12, "TC": 13, "TG": 14, "TT": 15}


def get_interest_values(cnx, exon_list, target_column, nt):
    """
    Return the list of ``target_column`` values.

    :param cnx: (sqlite3 connector object) connection to sed database
    :param exon_list: (list of 2 int) list of exons identified by the gene_id and exon position
    :param target_column: (string) the feature of interes
    :param nt: (string or NoneType) the nucleotide of interest
    :return: (list of float) the list of values of interest
    """
    if nt:
        if "mean_intron" in target_column:
            valuesa = np.array(
                figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, "iupac_upstream_intron", nt),
                dtype=float)
            valuesb = np.array(
                figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, "iupac_downstream_intron", nt),
                dtype=float)
            values = np.array([np.nanmedian([valuesa[i], valuesb[i]]) for i in range(len(valuesa))],
                              dtype=float)
        elif "introns" in target_column:
            valuesa = np.array(
                figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, "iupac_upstream_intron", nt),
                dtype=float)
            valuesb = np.array(
                figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, "iupac_downstream_intron", nt),
                dtype=float)
            values = np.concatenate((valuesa, valuesb))
        else:
            values = np.array(figure_producer.get_redundant_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt),
                              dtype=float)
    else:
        if target_column in ["median_flanking_intron_size", "min_flanking_intron_size", "introns_size"]:
            values_up = np.array(
                figure_producer.get_redundant_list_of_value(cnx, exon_list, "upstream_intron_size"),
                dtype=float)
            values_down = np.array(
                figure_producer.get_redundant_list_of_value(cnx, exon_list, "downstream_intron_size"),
                dtype=float)
            if target_column == "median_flanking_intron_size":
                values = np.array([np.nanmedian([values_up[i], values_down[i]]) for i in range(len(values_up))],
                                  dtype=float)
            elif target_column == "min_flanking_intron_size":
                values = np.array([np.nanmin([values_up[i], values_down[i]]) for i in range(len(values_up))],
                                  dtype=float)
            else:
                values = np.concatenate((values_up, values_down))
        else:
            values = np.array(figure_producer.get_redundant_list_of_value(cnx, exon_list, target_column),
                              dtype=float)
    return values


def get_control_exon(cnx, exon_type):
    """
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return: (list of 2 int) the list of control exons considerated
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos
                   FROM sed
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    print("control exon considerated : %s" % len(result))
    return result


def get_relative_value_of_a_project_or_sf(cnx, exon_list, target_column, control_dic, nt, operation, representation):
    """
    Return the median value of ``target_column`` in ``regulation`` exons for one exon list
    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param target_column:  (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param exon_list: (list of 2 int) list of exons (gene_id + exon_position)
    :param nt:  (string) the nt of interest
    :param control_dic: (dictionnary of list of float) median value of each possible control exons of \
    each feature in sed database.
    :param operation : (string) mean or median
    :param representation: (string) relative or aboslute
    :return: (float) the relative median value of the ``target_column`` in the list of exon comapred to control
    """
    values = get_interest_values(cnx, exon_list, target_column, nt)
    median_obs = eval("np.%s(values[~np.isnan(values)])" % operation)
    if representation == "absolute":
        return median_obs
    else:
        if nt:
            final_value = float(median_obs - control_dic[target_column][nt]) / control_dic[target_column][nt] * 100
        else:
            final_value = float(median_obs - control_dic[target_column]) / control_dic[target_column] * 100
        return final_value


def get_median_value(cnx, id_projects_sf_name, target_column, control_dic, regulation, operation,
                     representation, nt=None):
    """
    Return the median value of target_column in ``regulation`` exons of every  ``id_projects_sf_name``.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param id_projects_sf_name:  (list of int) the splicing lore id projects of every projects of interest
    :param target_column: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param control_dic: (dictionnary of list of float) median value of each possible control exons of \
    each feature in sed database.
    :param regulation: (list of string) up or down or up + down
    :param operation: (string) mean or median
    :param representation: (string) relative or absolute
    :param nt: (string) the nt of interest
    :return: (float) the relative median value (compared to control exons)  of ``target_column`` for every \
    ``regulation`` exons in every projects ``id_projects``
    """
    values_list = []
    try:
        int(id_projects_sf_name[0])
        sf_type = "project"
    except ValueError:
        sf_type = "sf"
    if sf_type == "project":
        for i in range(len(id_projects_sf_name)):
            exon_list = figure_producer.get_ase_events(cnx, id_projects_sf_name[i], regulation)
            final_value = get_relative_value_of_a_project_or_sf(cnx, exon_list, target_column, control_dic, nt,
                                                                operation, representation)
            values_list.append(final_value)
    else:
        for sf_name in id_projects_sf_name:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
            final_value = get_relative_value_of_a_project_or_sf(cnx, exon_list, target_column, control_dic, nt,
                                                                operation, representation)
            values_list.append(final_value)

    return values_list


def get_exons_values(cnx, sf_list, target_column1, target_column2, regulation):
    """
    Return the values of target_column in every`\
    `regulation`` exons regulated by a splicing factor in (one or multiple) cell lines.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param target_column1: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param target_column2: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    each feature in sed database.
    :param regulation: (list of string) up or down or up + down
    :return: 3 lists :

        * values : (list of list of float) each sublist corresponds to the value of `` target_column`` for \
        every exons regulated by a splicing factor
        * exon_name : (list of list of string) each sublist corresponds to the name of \
        every exons regulated by a splicing factor - the value in the sublist **i** position **j** \
         in the ``value`` and ``exon_name`` corresponds to the same exons
        * all_sf (list of string) list of each sf studied
    """
    if "$" in target_column1:
        target_column1, nt1 = target_column1.split("$")
    else:
        nt1 = None
    if "$" in target_column2:
        target_column2, nt2 = target_column2.split("$")
    else:
        nt2 = None
    exon_list = []
    if isinstance(sf_list[0], str):
        for sf_name in sf_list:
            exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
        exon_list = union_dataset_function.washing_events_all(exon_list)
    else:
        exon_list = sf_list
    print(len(exon_list))
    exon_name = ["%s_%s" % (union_dataset_function.get_gene_name(cnx, a[0]), a[1]) for a in exon_list]
    values1 = get_interest_values(cnx, exon_list, target_column1, nt1)
    values2 = get_interest_values(cnx, exon_list, target_column2, nt2)
    if len(exon_name) * 2 == len(values1):
        exon_name = ["%s_upstream" % a for a in exon_name ] + ["%s_downstream" % a for a in exon_name ]
    return values1, values2, exon_name


def get_gene_values(cnx, sf_list, target_column1, target_column2, regulation):
    """
    Return the values of target_column in every`\
    `regulation`` exons regulated by a splicing factor in (one or multiple) cell lines.

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param target_column1: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    :param target_column2: (string) the value for which we want to get the median value for the ``regulation`` \
    exon.
    each feature in sed database.
    :param regulation: (list of string) up or down or up + down
    :return: 3 lists :

        * values : (list of list of float) each sublist corresponds to the value of `` target_column`` for \
        every exons regulated by a splicing factor
        * exon_name : (list of list of string) each sublist corresponds to the name of \
        every exons regulated by a splicing factor - the value in the sublist **i** position **j** \
         in the ``value`` and ``exon_name`` corresponds to the same exons
        * all_sf (list of string) list of each sf studied
    """
    if "$" in target_column1:
        target_column1, nt1 = target_column1.split("$")
    else:
        nt1 = None
    if "$" in target_column2:
        target_column2, nt2 = target_column2.split("$")
    else:
        nt2 = None
    exon_list = []
    if isinstance(sf_list[0], str):
        for sf_name in sf_list:
            exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
        exon_list = union_dataset_function.washing_events_all(exon_list)
    else:
        exon_list = sf_list
    gene_id = []
    for val in exon_list:
        if val[0] not in gene_id:
            gene_id.append(val[0])
    gene_name = [union_dataset_function.get_gene_name(cnx, my_id) for my_id in gene_id]
    if nt1:
        values1 = get_list_of_value_iupac_dnt(cnx, exon_list, target_column1, nt1)
    else:
        values1 = figure_producer.get_list_of_value(cnx, exon_list, target_column1)

    if nt2:
        values2 = get_list_of_value_iupac_dnt(cnx, exon_list, target_column2, nt2)
    else:
        values2 = figure_producer.get_list_of_value(cnx, exon_list, target_column2)

    return values1, values2, gene_name


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
                   "SNRNP70", "SNRPC", "U2AF1", "U2AF2"]
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


def figure_creator(values_xaxis, values_yaxis, projects_names, regulation, name_xaxis, name_yaxis, ctrl, output,
                   name_fig, representation):
    """
    Create a scatter plot showing the potential correlation between projects.

    :param values_xaxis: (list of float) the list of relative median value displayed in xaxis
    :param values_yaxis: (list of float) the list of relative median value displayed in yaxis
    :param projects_names: (list of string) list of every projects name studied.
    :param regulation: (string) the regulation chosen
    :param name_xaxis: (string) the name of the feature represented in the x axis
    :param name_yaxis: (string) the name of the feature represented in the y axis
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    :param name_fig: (string) the name of the figure
    :param representation: (string) relative or absolute
    """
    nt1 = ""
    nt2 = ""
    if "$" in name_xaxis:
        name_xaxis, nt1 = name_xaxis.split("$")
        name_xaxis = name_xaxis.replace("iupac", "%s_nt" % nt1)
    if "$" in name_yaxis:
        name_yaxis, nt2 = name_yaxis.split("$")
        name_yaxis = name_yaxis.replace("iupac", "%s_nt" % nt2)
    if nt1 != "":
        x_name = name_xaxis.replace("iupac", "%s_nt" % nt1)
    else:
        x_name = name_xaxis
    if nt2 != "":
        y_name = name_yaxis.replace("iupac", "%s_nt" % nt1)
    else:
        y_name = name_yaxis
    if nt1 != nt2 and nt1 != "" and nt2 != "":
        print("warning : nt on x and y axis are different")
    if nt1 != "":
        nt = nt1
    else:
        nt = nt2
    # trace_pattern = color_maker(projects_names, values_xaxis, values_yaxis)
    data = []
    slope, intercept, r_value, p_value, std_err = stats.linregress(values_xaxis, values_yaxis)
    line = slope * np.array(values_xaxis) + intercept
    p2 = go.Scatter(x=values_xaxis,
                    y=line,
                    mode='lines',
                    line=dict(color='red', width=3),
                    name="Fit"
                    )
    data.append(p2)
    # for key in trace_pattern:
    #     data.append(go.Scatter(
    #         x=trace_pattern[key][0],
    #         y=trace_pattern[key][1],
    #         name=key,
    #         mode='markers',
    #         text=trace_pattern[key][3],
    #         marker=dict(
    #             size=10,
    #             color=trace_pattern[key][2][0],
    #             line=dict(width=1))
    #     ))
    data.append(go.Scatter(
            x=values_xaxis,
            y=values_yaxis,
            name="projects",
            mode='markers',
            text=projects_names,
            marker=dict(
                size=10,
                line=dict(width=1))))
    cor, pval = stats.pearsonr(values_xaxis, values_yaxis)

    main_title = 'Correlation between %s and %s frequency in genes containing %s-regulated exons ' \
                 'in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' \
                 % (x_name, y_name, regulation, ctrl, round(cor, 2), pval)
    x_title = '%s %s' % (representation, x_name)
    y_title = '%s %s frequency' % (representation, y_name)
    figname = '%s%s_%s_%s.html' % (output, nt, name_fig, regulation)

    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title.replace("absolute", ""),
            tickfont=dict(size=25)),
        yaxis=dict(
            title=y_title.replace("absolute", ""),
            tickfont=dict(
                size=25)),
        showlegend=True
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname,
                        auto_open=False)


def get_relative_values(list_values, ctrl_dic, name_list, representation):
    """
    Get the relative value of a list given.

    :param list_values: (list of float) the list of value
    :param name_list: (string) the name of the list
    :param ctrl_dic: (dictionary) median for every exon within ``ctrl_dic``
    :param representation: (string) relative or absolute
    :return: - (list of float) the list of relative values toward control exons
             - (string) the new name for the list
    """
    if max(list_values) > 1000:
        if representation == "relative":
            if "$" in name_list:
                target_column, nt = name_list.split("$")
                list_values2 = [float(math.log10(val) - math.log10(ctrl_dic[target_column][nt])) /
                                math.log10(ctrl_dic[target_column][nt]) * 100 for val in list_values]
            else:
                # for i in range(len(exon_name)):
                #     if exon_name[i] == "CDC25C_6_downstream":
                #         print(exon_name[i], name_list)
                #         print("float(%s - %s) / %s  * 100 " % (list_values[i], ctrl_dic[name_list], ctrl_dic[name_list]))
                list_values2 = [float(math.log10(val) - math.log10(ctrl_dic[name_list])) /
                                math.log10(ctrl_dic[name_list]) * 100 for val in list_values]

        else:
            list_values2 = [math.log10(val) for val in list_values]
        name_list = "log_%s" % name_list
    else:
        if representation == "relative":
            if "$" in name_list:
                target_column, nt = name_list.split("$")
                list_values2 = [float(val - ctrl_dic[target_column][nt]) /
                                ctrl_dic[target_column][nt] * 100 for val in list_values]
            else:
                list_values2 = [float(val - ctrl_dic[name_list]) / ctrl_dic[name_list] * 100 for val in list_values]
        else:
            list_values2 = list_values
    return list_values2, name_list


def figure_creator_exon(values_xaxis, values_yaxis, regulation, name_xaxis, name_yaxis, ctrl,
                        exon_name, output, name_fig, representation):
    """
    Create a scatter plot showing the potential correlation for every exons regulated each splicing factor in \
    multiple cell lines.

    :param values_xaxis: (list of float) the list of relative value of a particular feature for \
    exons ``regulation`` regulated by ``sf_name``
    :param values_yaxis: (list of float) the list of relative value of the frequency in ``nt`` for \
    exons ``regulation`` regulated by ``sf_name``
    :param regulation: (list of string) the regulation chosen up down or up and down
    :param name_xaxis: (string) the name of the features displayed in xaxis
    :param name_yaxis: (string) the name of the features displayed in yaxis
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    :param exon_name: (list of string) the name of the interest exons
    :param name_fig: (string) the name of the graphics
    :param representation: (string) relative or absolute
    """
    data = []
    nt1 = ""
    nt2 = ""
    if "$" in name_xaxis:
        name_xaxis, nt1 = name_xaxis.split("$")
        name_xaxis = name_xaxis.replace("iupac", "%s_nt" % nt1)
    if "$" in name_yaxis:
        name_yaxis, nt2 = name_yaxis.split("$")
        name_yaxis = name_yaxis.replace("iupac", "%s_nt" % nt2)
    if nt1 != "":
        x_name = name_xaxis.replace("iupac", "%s_nt" % nt1)
    else:
        x_name = name_xaxis
    if nt2 != "":
        y_name = name_yaxis.replace("iupac", "%s_nt" % nt1)
    else:
        y_name = name_yaxis
    if nt1 != nt2 and nt1 != "" and nt2 != "":
        print(nt1, nt2)
        print("warning : nt on x and y axis are different")
    if nt1 != "":
        nt = nt1
    else:
        nt = nt2
    slope, intercept, r_value, p_value, std_err = stats.linregress(values_xaxis, values_yaxis)
    line = slope * np.array(values_xaxis) + intercept
    cor, pval = stats.pearsonr(values_xaxis, values_yaxis)
    data.append(go.Scatter(x=values_xaxis,
                           y=line,
                           mode='lines',
                           line=dict(width=3,
                                     color="red"),
                           name="Fit : p=%.2E, c=%.2f" % (pval, cor)
                           ))
    data.append(go.Scatter(
        x=values_xaxis,
        y=values_yaxis,
        name="exon down of every projects",
        mode='markers',
        text=exon_name,
        marker=dict(
            size=7,
            color="navy"
            )
    ))

    main_title = 'Correlation between %s and %s frequency in genes containing %s-regulated exons ' \
                 'in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' \
                 % (x_name, y_name, regulation, ctrl, round(cor, 2), pval)

    x_title = '%s %s' % (representation, x_name)
    y_title = '%s %s frequency' % (representation, y_name)
    figname = '%s%s_%s_%s_%s.html' % (output, nt, name_fig, regulation, representation)
    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title.replace("absolute", ""),
            tickfont=dict(size=25),
            showgrid=False),
        yaxis=dict(
            title=y_title.replace("absolute", ""),
            tickfont=dict(size=25),
            showgrid=False),
        showlegend=True
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname, auto_open=False)


def density_creator_exon(values_xaxis, values_yaxis, regulation, name_xaxis, name_yaxis, ctrl, output, name_fig,
                         representation):
    """
    Create a density plot showing the potential correlation for every exons regulated each splicing factor in \
    multiple cell lines.

    :param values_xaxis: (list of float) the list of relative value of a particular feature for \
    exons ``regulation`` regulated by ``sf_name``
    :param values_yaxis: (list of float) the list of relative value of the frequency in ``nt`` for \
    exons ``regulation`` regulated by ``sf_name``
    :param regulation: (list of string) the regulation chosen up down or up and down
    :param name_xaxis: (string) the name of the features displayed in xaxis
    :param name_yaxis: (string) the name of the features displayed in yaxis
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    :param name_fig: (string) the name of the graphics
    :param representation: (string) relative or absolute
    """
    data = []
    nt1 = ""
    nt2 = ""
    if "$" in name_xaxis:
        name_xaxis, nt1 = name_xaxis.split("$")
        name_xaxis = name_xaxis.replace("iupac", "%s_nt" % nt1)
    if "$" in name_yaxis:
        name_yaxis, nt2 = name_yaxis.split("$")
        name_yaxis = name_yaxis.replace("iupac", "%s_nt" % nt2)
    if nt1 != "":
        x_name = name_xaxis.replace("iupac", "%s_nt" % nt1)
    else:
        x_name = name_xaxis
    if nt2 != "":
        y_name = name_yaxis.replace("iupac", "%s_nt" % nt1)
    else:
        y_name = name_yaxis
    if nt1 != nt2 and nt1 != "" and nt2 != "":
        print(nt1, nt2)
        print("warning : nt on x and y axis are different")
    if nt1 != "":
        nt = nt1
    else:
        nt = nt2
    slope, intercept, r_value, p_value, std_err = stats.linregress(values_xaxis, values_yaxis)
    rep = values_xaxis
    if representation == "relative":
        line = slope * np.array(values_xaxis) + intercept
    else:
        if intercept < 0:
            rep = [-intercept / slope] + list(values_xaxis)
            line = slope * np.array(rep) + intercept
        else:
            rep = [0] + list(values_xaxis)
            line = slope * np.array(rep) + intercept
    cor, pval = stats.pearsonr(values_xaxis, values_yaxis)
    data.append(go.Scatter(x=rep,
                           y=line,
                           mode='lines',
                           line=dict(width=3),
                           name="Fit : p=%.2E, c=%.2f" % (pval, cor),
                           showlegend=False
                           ))
    data.append(dict(type='histogram2dcontour',
                     x=values_xaxis, y=values_yaxis, ncontours=50,
                     # xbins=dict(start=min(values_xaxis)-bins, end=max(values_xaxis)+bins, size=bins),
                     # ybins=dict(start=min(values_yaxis)-bins, end=max(values_yaxis)+bins, size=bins),
                     colorscale='Hot', reversescale=True, contours={'showlines': False}
                     )
                )

    main_title = 'Density between %s and %s frequency in genes containing %s-regulated exons ' \
                 'in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' \
                 % (x_name, y_name, regulation, ctrl, round(cor, 2), pval)
    x_title = '%s %s' % (representation, x_name)
    y_title = '%s %s frequency' % (representation, y_name)
    figname = '%sdensity_%s_%s_%s_%s.html' % (output, nt, name_fig, regulation, representation)
    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title.replace("absolute", ""),
            showgrid=False),
        yaxis=dict(
            title=y_title.replace("absolute", ""),
            showgrid=False),
        showlegend=True,
        shapes=[dict(type='line', x0=0, y0=min(min(values_yaxis), 0), x1=0, y1=max(values_yaxis) * 1.4, line=dict(color="black")),
                dict(type='line', x0=min(min(values_xaxis), 0), y0=0, x1=max(values_xaxis) * 1.4, y1=0, line=dict(color="black"))]
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname, auto_open=False)


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
    values_list1 = list(values_list1)
    indice_2_remove = []
    list_exon_name = np.array(list_exon_name)
    for i in range(len(values_list1)):
        if values_list1[i] is None or values_list2[i] is None or np.isnan(values_list1[i]) or np.isnan(values_list2[i]):
            indice_2_remove.append(i)
    values_list1 = np.delete(values_list1, indice_2_remove)
    values_list2 = np.delete(values_list2, indice_2_remove)
    list_exon_name = np.delete(list_exon_name, indice_2_remove)
    return values_list1, values_list2, list(list_exon_name)


def define_couple_targets(xaxis, yaxis, nt_list):
    """
    Create a list of list of 2 element, each sublist will indicate the unit that will be represent in the axis \
    of the correlation graphics

    :param xaxis: (string) the value of xaxis
    :param yaxis: (string) the value of yaxis
    :param nt_list: (string) a list of nucleotides
    :return: (list of list of 2 strings)  each sublist corresponds to the value in the x and y axis \
    that will be represented in the graphics
    """
    couple_targets = []
    if "iupac" not in xaxis and "iupac" not in yaxis:
        couple_targets.append([xaxis, yaxis])
        return couple_targets
    if "iupac" in xaxis and "iupac" not in yaxis:
        for nt in nt_list:
            couple_targets.append(["%s$%s" % (xaxis, nt), yaxis])
        return couple_targets
    if "iupac" not in xaxis and "iupac" in yaxis:
        for nt in nt_list:
            couple_targets.append([xaxis, "%s$%s" % (yaxis, nt)])
        return couple_targets
    if "iupac" in xaxis and "iupac" in yaxis:
        for nt in nt_list:
            couple_targets.append(["%s$%s" % (xaxis, nt), "%s$%s" % (yaxis, nt)])
        return couple_targets


def get_axis_value(target, cnx, id_projects_sf_name, ctrl_dic, regulation, operation, representation):
    """
    Get the wanted list of value

    :param target: (string) a iupac of non iupac target
    :param cnx: (sqlite3 connect object) connection to seddb
    :param id_projects_sf_name: (list of int or string) list of every splicing lore project that we want to represent \
    as point in our figure or list of every sf name that we want to analyze
    :param ctrl_dic: (dictionary) median for every exon within ``ctrl_dic``
    :param regulation: (string) down or up
    :param operation: (string) median or mean
    :param representation: (string) relative or absolute
    :return: (the list of float) list of median value of the feature ``target`` for each exon regulated in each project
    """
    if "$" not in target:
        value_list = get_median_value(cnx, id_projects_sf_name, target, ctrl_dic, regulation, operation,
                                      representation, nt=None)
    else:
        target, nt = target.split("$")
        value_list = get_median_value(cnx, id_projects_sf_name, target, ctrl_dic, regulation, operation,
                                      representation, nt=nt)
    return value_list


def get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt_dnt: (string) a nucleotide or di_nucleotide
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                if len(nt_dnt) == 1:
                    res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                else:
                    res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
            else:
                res.append(None)
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                if r is not None:
                    if len(nt_dnt) == 1:
                        res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                    else:
                        res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
                else:
                    res.append(None)
                redundancy_gene_dic[exon[0]] = 1
    return res


def remove_wrong_size_values(list_values1, name_value1, name_exon, list_values2, name_value2):
    """
    Remove the size value of 0 and below it.

    :param list_values1: (list of float) list of values for the feature ``name_values1``
    :param name_value1: (string) the name of the feature1
    :param name_exon: (list of string) the name of each exons for which we have the the values in list_values1 and 2
    :param list_values2:(list of float) list of values for the feature ``name_values2``
    :param name_value2: (string) the name of the feature2
    :return: the list of values without 0 and negative values
    """
    if "size" not in name_value2 and "size" not in name_value1:
        return list_values1, list_values2, name_exon
    else:
        if "size" in name_value1:
            i = 0
            while i < len(list_values1):
                if list_values1[i] < 1:
                    del(list_values1[i])
                    del(list_values2[i])
                    del(name_exon[i])
                    i -= 1
                i += 1
        if "size" in name_value2:
            i = 0
            while i < len(list_values2):
                if list_values2[i] < 1:
                    del(list_values1[i])
                    del(list_values2[i])
                    del(name_exon[i])
                    i -= 1
                i += 1
        return list_values1, list_values2, name_exon


def project_t2_remove(cnx, sf_list):
    """
    Get the id of the projects we want to remove
    :param cnx: (sqlite3 connect object) connection to sed database
    :param sf_list: (string) list of splicing factor
    :return: (list of int) the id project to remove
    """
    cursor = cnx.cursor()
    list_projects = []
    query = "SELECT id FROM rnaseq_projects WHERE sf_name = '%s'"
    for sf in sf_list:
        cursor.execute(query % sf)
        res = cursor.fetchall()
        for project in res:
            list_projects.append(project[0])
    return list_projects


def removing_projects(full_projects_list, name_projects_list, project_2_remove):
    """
    Return ``full_projects_list`` without every project present in ``project_2_remove``

    :param full_projects_list: (list of int) list of every project in splicing lore
    :param name_projects_list:  (list of string) the nam of every project in splicing lore
    :param project_2_remove: (list of int) list of every project we need to remove in splicing lore
    :return: (2 lists):
        * (list of int) the list of project we want to analye
        * (list of string) the list of projects name we want to  analyse
    """
    wanted_project_id = []
    wanted_project_name = []
    for i in range(len(full_projects_list)):
        if full_projects_list[i] not in project_2_remove:
            wanted_project_id.append(full_projects_list[i])
            wanted_project_name.append(name_projects_list[i])
    return wanted_project_id, wanted_project_name


def difference(cnx, list1, list2, regulation):
    """
    Return the exons regulated by the factors in list1 if they are not regulated by the factors in list2
    :param cnx: (sqlite3 connect object) connection to sed database
    :param list1: (list of string) list of splicing factors
    :param list2:  (list of strings) list of splicing factros
    :param regulation: (string) the exons with the regulation ``regulation`` regulated by the splicing factors in \
    ``list1`` or ``list2``
    :return:(list of list of 2 int
    """
    exon_list1 = []
    exon_list2 = []
    for sf_name in list1:
        exon_list1 += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_list1 = union_dataset_function.washing_events_all(exon_list1)
    for sf_name in list2:
        exon_list2 += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_list2 = union_dataset_function.washing_events_all(exon_list2)
    return [exon for exon in exon_list1 if exon not in exon_list2]


def get_concidered_exons(cnx, exon_class, regulation):
    """
    Get the exon considered
    :param cnx: (sqlite3 connect object) connection to sed databse
    :param exon_class: (string) None, GC, AT or GC-AT, the exons that we consider, if \
    exon_class is set to None, all exon regulated by a splicing factor are considered
    :param regulation: (string) up or down
    :return: (list of string/list 2 int) list of sf_name if exon_class is None or list of exons if \
                     exon_class is not None.
    """
    if exon_class is None:
        list_sf = group_factor.get_wanted_sf_name(None)
        return list_sf
    else:
        if exon_class == "AT":
            exon_list = difference(cnx, group_factor.at_rich_down, group_factor.gc_rich_down, regulation)
        elif exon_class == "GC":
            exon_list = difference(cnx, group_factor.gc_rich_down, group_factor.at_rich_down, regulation)
        elif exon_class == "GC-AT":
            exon_list = difference(cnx, group_factor.at_rich_down, group_factor.gc_rich_down, regulation)
            exon_list += difference(cnx, group_factor.gc_rich_down, group_factor.at_rich_down, regulation)
            print("exons GC-AT concidered : %s" % len(exon_list))
        else:
            exon_list = get_control_exon(cnx, exon_class)
        return exon_list



def main(level, xaxis, yaxis, name_fig, exon_type, nt_list, exon_class, operation, representation):
    """
    Create the correlation matrix (gene_size vs iupac)
    """
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic, full_dic = control_exon_adapter.control_handler(cnx, exon_type, operation)
    nt_list = nt_list.split(",")
    regulation = "down"
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/new_correlation/"
    couple_targets = define_couple_targets(xaxis, yaxis, nt_list)
    if level == "project":
        name_fig += "_project"
        id_projects, name_projects = group_factor.get_id_and_name_project_wanted(cnx, None)
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        for i in range(len(couple_targets)):
            values_xaxis = get_axis_value(couple_targets[i][0], cnx, id_projects, ctrl_dic, regulation, operation,
                                          representation)
            values_yaxis = get_axis_value(couple_targets[i][1], cnx, id_projects, ctrl_dic, regulation, operation,
                                          representation)
            if len(values_xaxis) != len(values_yaxis):
                print("Warning the list of value don't have the same length")
            figure_creator(values_xaxis, values_yaxis, name_projects, regulation, couple_targets[i][0],
                           couple_targets[i][1], exon_type, output,
                           name_fig, representation)
    elif level == "exon":
        if exon_class is None:
            list_sf = group_factor.get_wanted_sf_name(None)
        else:
            exon_list = get_concidered_exons(cnx, exon_class, regulation)
            name_fig += "_%s_exons" % exon_class
        name_fig += "_Exon_LVL"
        if not os.path.isdir(output):
            os.mkdir(output)
        gene_columns = ["median_intron_size", "iupac_gene", "gene_size"]
        if xaxis in  gene_columns and yaxis in gene_columns:
            for i in range(len(couple_targets)):
                if exon_class is None:
                    value_xaxis, value_yaxis, exon_name = get_gene_values(cnx, list_sf, couple_targets[i][0],
                                                                          couple_targets[i][1], regulation)
                else:
                    value_xaxis, value_yaxis, exon_name = get_gene_values(cnx, exon_list, couple_targets[i][0],
                                                                          couple_targets[i][1], regulation)
                value_xaxis, value_yaxis, exon_name =  remove_none_values(value_xaxis, value_yaxis, exon_name)
                value_xaxis, name_x = get_relative_values(value_xaxis, ctrl_dic, couple_targets[i][0], representation)
                value_yaxis, name_y = get_relative_values(value_yaxis, ctrl_dic, couple_targets[i][1], representation)
                figure_creator_exon(value_xaxis, value_yaxis, regulation, name_x, name_y, exon_type,
                                    exon_name, output, name_fig, representation)
                density_creator_exon(value_xaxis, value_yaxis, regulation, name_x, name_y, exon_type, output, name_fig,
                                     representation)
        else:
            for i in range(len(couple_targets)):
                if exon_class is None:
                    value_xaxis, value_yaxis, exon_name = get_exons_values(cnx, list_sf, couple_targets[i][0],
                                                                           couple_targets[i][1], regulation)
                else:
                    value_xaxis, value_yaxis, exon_name = get_exons_values(cnx, exon_list, couple_targets[i][0],
                                                                           couple_targets[i][1], regulation)
                value_xaxis, value_yaxis, exon_name = remove_none_values(value_xaxis, value_yaxis, exon_name)
                value_xaxis, value_yaxis, exon_name = remove_wrong_size_values(value_xaxis, couple_targets[i][0],
                                                                               exon_name, value_yaxis,
                                                                               couple_targets[i][1])

                if len(value_xaxis) != len(value_yaxis) or len(value_xaxis) != len(exon_name):
                    print("Warning the list of value do'nt have the same length")
                value_xaxis, name_x = get_relative_values(value_xaxis, ctrl_dic, couple_targets[i][0], representation)
                value_yaxis, name_y = get_relative_values(value_yaxis, ctrl_dic, couple_targets[i][1], representation)
                #if exon_class in ["GC", "AT", "GC-AT"]:
                figure_creator_exon(value_xaxis, value_yaxis, regulation, name_x, name_y, exon_type,
                                        exon_name, output, name_fig, representation)
                density_creator_exon(value_xaxis, value_yaxis, regulation, name_x, name_y, exon_type, output, name_fig,
                                     representation)
    else:
        name_fig += "_union"
        list_sf = group_factor.get_wanted_sf_name(None)
        if not os.path.isdir(output):
            os.mkdir(output)
        for i in range(len(couple_targets)):
            values_xaxis = get_axis_value(couple_targets[i][0], cnx, list_sf, ctrl_dic, regulation, operation,
                                          representation)
            values_yaxis = get_axis_value(couple_targets[i][1], cnx, list_sf, ctrl_dic, regulation, operation,
                                          representation)
            if len(values_xaxis) != len(values_yaxis):
                print("Warning the list of value don't have the same length")
            figure_creator(values_xaxis, values_yaxis, list_sf, regulation, couple_targets[i][0],
                           couple_targets[i][1], exon_type, output,
                           name_fig, representation)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Create a correlation graphics

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("required argument")

    parser.add_argument('--level', dest='level',
                        help="""exon for each point to be an exon, 
                                project for each point tb be a project and sf for each point to be a sf""",
                        default="project")

    parser.add_argument('--name', dest='name', help='the name of the graphics', default="correlation")

    parser.add_argument('--nt_list', dest='nt_list', help='the list of nucleotide you '
                                                          'want to analyse (they must be separated by a coma)',
                        default="S,W")

    parser.add_argument('--exon_type', dest='exon_type', help='name of control exons', default="CCE")
    parser.add_argument('--exon_class', dest='exon_class', help='the class of exon we want to study', default=None)
    parser.add_argument("--operation", dest="operation", help="the type of heatmap you want to produce (mean/median)",
                        default="median")
    parser.add_argument("--representation", dest="representation", help="relative to represent the relative frequencies"
                                                                        "of exons of absolute else", default="relative")

    required_args.add_argument('--xaxis', dest='xaxis',
                               help="""element of the xaxis""",
                               required=True)

    required_args.add_argument('--yaxis', dest='yaxis',
                               help="""element in the yaxis""",
                               required=True)
    args = parser.parse_args()  # parsing arguments
    # executing the program
    class_approuved = ["AT", "GC", "GC-AT", "CCE", "ACE", "ALL", None]
    if args.exon_class not in class_approuved:
        parser.error("ERROR : wrong value for --exon_class argument : only %s are allowed" % " ".join(class_approuved))
    if args.representation not in ["absolute", "relative"]:
        parser.error("Error : the argument representation can only by either absolute or relative")
    main(args.level, args.xaxis, args.yaxis, args.name, args.exon_type, args.nt_list, args.exon_class, args.operation,
         args.representation)


if __name__ == "__main__":
    launcher()
