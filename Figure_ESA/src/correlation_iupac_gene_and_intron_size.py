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


def figure_creator(values_size, values_iupac, projects_names, regulation, size_scale, nt_name, ctrl, output):
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
    layout = go.Layout(
        title='Correlation between %s and %s frequency for %s exons in every splicing lore project<br>'
              '(relative value against %s control) - cor : %s - pval : %.2E'
              % (size_scale, nt_name, regulation, ctrl, round(cor, 2), pval),
        hovermode='closest',
        xaxis=dict(
            title='relative %s' % size_scale),
        yaxis=dict(
            title='relative %s frequency' % nt_name),
        showlegend=True
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename='%s%s_%s_correlation_graphs_%s.html'
                                      % (output, nt_name, size_scale, regulation[0]),
                        auto_open=False)


def main():
    """
    Create the correlation matrix.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    id_projects, name_projects = figure_producer.get_interest_project(cnx)
    nt_list = ["A", "C", "G", "T", "S", "W", "Y", "R"]
    output_niv0 = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/correlation_size_frequency/"
    # If the output directory does not exist, then we create it !
    if not os.path.isdir(output_niv0):
        os.mkdir(output_niv0)
    regulations = [["up"], ["down"], ["up-down"]]
    targets1 = ["gene_size", "median_intron_size"]
    target2 = "iupac_gene"
    for my_regulation in regulations:
        output_niv1 = output_niv0 + my_regulation[0] + "/"
        if not os.path.isdir(output_niv1):
            os.mkdir(output_niv1)
        for target1 in targets1:
            value_target1 = get_median_value(cnx, id_projects, target1, ctrl_dic, my_regulation, nt=None)
            for nt in nt_list:
                value_target2 = get_median_value(cnx, id_projects, target2, ctrl_dic, my_regulation, nt=nt)
                figure_creator(value_target1, value_target2, name_projects, my_regulation[0],
                               target1, nt, exon_type, output_niv1)


if __name__ == "__main__":
    main()
