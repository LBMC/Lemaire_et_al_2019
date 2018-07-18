#!/usr/bin/env python3

# -*- coding: utf-8

# Import modules
import get_intron_information
import sqlite3
import sys
import os
import numpy as np
import figure_producer
import plotly
import plotly.graph_objs as go
from scipy import stats


# function
def database_connection(db_name):
    """
    :param db_name: (string) the name of the database to connect

    Allow connection to the database ``db_name``
    """
    return sqlite3.connect(db_name)


def exon_finder(cnx_fasterdb, exon_type):
    """
    Find for every exons in fasterDB lite, its official gene symbol, its \
    gene id and its position within a gene.

    :param cnx_fasterdb: (sqlite3 object) allows connection to fasterDB lite
    :param exon_type: the type of control exon we want to use.
    :return: (list of tuple of string int int). Each sublist of the list \
    contains:
        * the gene symbol of the gene that contains the exon
        * the gene id of the gene that contains the exon
        * the position of tne exon within the gene
    """
    cursor = cnx_fasterdb.cursor()
    if exon_type != "ALL":
        query = """SELECT t2.official_symbol , t1.id_gene, t1.pos_on_gene
                   FROM genes t2, exons t1
                   WHERE t2.id = t1.id_gene
                   AND exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT t2.official_symbol , t1.id_gene, t1.pos_on_gene
                   FROM genes t2, exons t1
                   WHERE t2.id = t1.id_gene"""
    cursor.execute(query)
    return cursor.fetchall()


def get_exon_info(cnx, info_list, debug, display=True):
    """
    Get every information we need on an exon and it's surrounding introns

    :param debug: (int) 0 no debug, 1 debug mode
    :param cnx: (sqlite3 object) return all the information we need to connect to FasterDB lite
    :param info_list: (list of list of string and int and int) each sublist contains \
    a string : gene_symbol and 2 int : the gene_id and the exon position on gene respectively
    :param display: (boolean) True if we want to display the progression of the process used to get exon's data.
    :return: (a list of ExonClass object) every object contains the frequencies of every nucleotide \
    in a exon sequences (exon in ``info_list``) and it's surrounding introns.
    """
    if display:
        print("Getting exons information !")
    exon_list = []
    get_intron_information.set_debug(debug)
    count = 0
    ll = str(len(info_list))
    for exon_info in info_list:
        exon_list.append(get_intron_information.ExonClass(cnx, exon_info[0], exon_info[1], exon_info[2]))
        count += 1
        percent = round(float(count) / len(info_list) * 100, 1)
        if display:
            sys.stdout.write("Progression : " + str(count) + " / " + ll + " - " + str(percent) + " %\r")
            sys.stdout.flush()
    return exon_list


def get_ase_events(cnx, id_project, regulation):
    """
    Get every exon up or down regulated in a particular project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :param regulation: (string)) up or down
    :return: (list of tuple of 1 string) 2 int) each sublist corresponds to an exon (gene_symbol, gene_id + \
    exon_position on gene)
    """
    if regulation == "up":
        regulation = ">= 0.1"
    else:
        regulation = "<= -0.1"
    cursor = cnx.cursor()
    query = """SELECT gene_symbol, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue_glm_cor <= 0.05""" % (id_project, regulation)
    cursor.execute(query)
    res = cursor.fetchall()
    if len(res) == 0:
            query = """SELECT gene_symbol, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue <= 0.05""" % (id_project, regulation)
            cursor.execute(query)
            res = cursor.fetchall()
    return res


def get_data_dictionary(exon_list):
    """
    Get the median values of each nucleotides (iupac) for every exons in ``exon_list`` \
    for exons, downstream and upstream intron sequences.

    :param exon_list: (list of ExonClass object) iupac dictionaries for exon, upstream and downstream introns.
    :return: (dictionary of 3 dictionaries of float) each dictionaries corresponds to the frequencies of \
    every iupac nucleotides (key) in :
     * exon sequences
     * upstream intron sequence
     * downstream intron sequences
    """
    nt_list = ["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M"]
    up_intron_iupac = {a: [] for a in nt_list}
    down_intron_iupac = {a: [] for a in nt_list}
    exon_iupac = {a: [] for a in nt_list}
    res_dic = {"exon_iupac": exon_iupac, "upstream_intron_iupac": up_intron_iupac,
               "downstream_intron_iupac": down_intron_iupac}
    for exon in exon_list:
        for nt in nt_list:
            if exon.iupac_exon is not None:
                res_dic["exon_iupac"][nt].append(exon.iupac_exon[nt])
            else:
                res_dic["exon_iupac"][nt].append(None)
            if exon.iupac_up_intron is not None:
                res_dic["upstream_intron_iupac"][nt].append(exon.iupac_up_intron[nt])
            else:
                res_dic["upstream_intron_iupac"][nt].append(None)
            if exon.iupac_down_intron is not None:
                res_dic["downstream_intron_iupac"][nt].append(exon.iupac_down_intron[nt])
            else:
                res_dic["downstream_intron_iupac"][nt].append(None)
    get_intron_information.printd("res")
    get_intron_information.printd(res_dic)
    for key in res_dic:
        for nt in nt_list:
            val = np.array(res_dic[key][nt], dtype=float)
            res_dic[key][nt] = round(np.median(val[~np.isnan(val)]), 1)
    get_intron_information.printd("median_res")
    get_intron_information.printd(res_dic)
    return res_dic


def get_project_data_list(cnx_sed, cnx, id_projects, regulation, debug):
    """
    Get the dictionary that contains the median frequency of every iupac nucleotides in every splicing lore \
    project for up or down regulated exons.

    :param cnx_sed: (sqlite3 connect object) allow connexion to sed database
    :param cnx: (sqlite3 connect object) allow connexion to fasterDB database
    :param id_projects: (lits of int) list of splicing lore id project
    :param regulation: (string) up if we want to work on  up-regulated exon, down if we want  to work on \
    down-regulated exons.
    :param debug: (int) 1 debug mode enabled 0 debug mode disabled
    :return: (list of dictionaries of 3 dicionnaries  of float) each dictionaries contains 3 keys (exon_iupac, \
    upstream_intron_iupac and downstream_intron_iupac) that corresponds to the region they corresponds to. \
    Each one of these keys are linked to another dictionary containing the frequency \
    of every iupac nucleotide in that region.
    """
    res = []
    for id_project in id_projects:
        exon_list = get_ase_events(cnx_sed, id_project, regulation)
        exon_data = get_exon_info(cnx, exon_list, debug, False)
        exon_dic = get_data_dictionary(exon_data)
        res.append(exon_dic)
    return res


def get_list_of_values(exon_list_dic, target_column, control_dic, nt):
    """
    Return a list of interest float values.

    :param exon_list_dic: (list of dictionaries of 3 dicionnaries  of float) each dictionaries contains 3 keys \
    (exon_iupac, upstream_intron_iupac and downstream_intron_iupac) that corresponds to the region they corresponds to.\
     Each one of these keys are linked to another dictionary containing the frequency \
    of every iupac nucleotide in thta region.
    :param target_column: (string) the region we want to work on (exon, upstream_intron, downstream_intron)
    :param control_dic: (dictionaries containing 3 dictionaries of float) the dictionary has 3 keys (exon_iupac,
    upstream_intron_iupac, downstream_intron_iupac). Each key is linked to a dictionary that links each iupac \
    nucleotide to it's frequency in the region corresponding to that key for very control exons.
    :param nt: (string) a nucleotide (A, C, G, T, S, W, R, Y, K, M)
    :return: (list of float) the list of relative median iupac frequencies (toward control) in the ``target_column`` \
    regions of every group of regulated exon in every splicing lore project.
    """
    list_value = []
    for dic_iupac in exon_list_dic:
        list_value.append((dic_iupac[target_column][nt] -
                           control_dic[target_column][nt]) /
                          control_dic[target_column][nt] * 100)
    return list_value


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


def figure_creator(values_abscissa, values_ordinate, projects_names, regulation, name_abscissa, name_ordinate,
                   nt_name, ctrl, output):
    """
    Create a scatter plot showing the potential correlation between projects.

    :param values_abscissa: (list of float) the list of relative median value of ``name_absicssa`` for \
    every projects named ``projects_names``.
    :param values_ordinate: (list of float) the list of relative median value of ``name_ordinate`` for \
    every projects named ``projects_names``.
    :param projects_names: (list of string) list of every projects name studied.
    :param regulation: (list of string) the regulation chosen up or down
    :param name_abscissa: (string) exon_iupac
    :param name_ordinate: (string) upstream_intron_iupac or downstream_intron_iupac
    :param nt_name: (string) the name of the nucleotide of interest
    :param ctrl: (string) the control exon used to calculate relative frequency.
    :param output: (string) path where the results will be created
    """
    trace_pattern = color_maker(projects_names, values_abscissa, values_ordinate)
    data = []
    slope, intercept, r_value, p_value, std_err = stats.linregress(values_abscissa, values_ordinate)
    line = slope * np.array(values_abscissa) + intercept
    p2 = go.Scatter(x=values_abscissa,
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
    cor, pval = stats.pearsonr(values_abscissa, values_ordinate)
    name_abscissa = name_abscissa.replace("_iupac", "")
    name_ordinate = name_ordinate.replace("_iupac", "")
    main_title = 'Correlation between %s and %s for %s frequency in %s exons in every splicing lore project<br> (relative value against %s control) - cor : %s - pval : %.2E' % (name_abscissa, name_ordinate, nt_name, regulation, ctrl, round(cor, 2), pval)
    x_title = 'relative %s frequency in %s' % (nt_name, name_abscissa)
    y_title = 'relative %s frequency in %s' % (nt_name, name_ordinate)
    figname = '%s%s_%s_%s_correlation_graphs_%s.html' % (output, nt_name, name_abscissa, name_ordinate, regulation)

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


def main():
    """
    Create correlation figure for iupac downstream and upstream intron vs iupac exon
    """
    ctrl = "CCE"
    nt_list = ["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M"]
    debug = 0  # 1 enabled, 0 disabled
    fasterdb = os.path.realpath(os.path.dirname("__file__")) + "/data/fasterDB_lite.db"
    seddb = os.path.realpath(os.path.dirname("__file__")) + "/data/sed.db"
    output = os.path.realpath(os.path.dirname("__file__")) + "/result/correlation_iupac_introns_exons/"
    if not os.path.isdir(output):
        os.mkdir(output)
    print(fasterdb)
    cnx = database_connection(fasterdb)
    cnx_sed = database_connection(seddb)
    exon_list = exon_finder(cnx, ctrl)
    exon_control = get_exon_info(cnx, exon_list, debug)
    control_dic = get_data_dictionary(exon_control)
    id_projects, name_projects = figure_producer.get_interest_project(cnx_sed)
    target1 = "exon_iupac"
    targets2 = ["upstream_intron_iupac", "downstream_intron_iupac"]
    for regulation in ["up", "down"]:
        exon_dic_list = get_project_data_list(cnx_sed, cnx, id_projects, regulation, debug)
        for target2 in targets2:
            final_output = output + target2 + "/"
            if not os.path.isdir(final_output):
                os.mkdir(final_output)
            for nt in nt_list:
                values_target1 = get_list_of_values(exon_dic_list, target1, control_dic, nt)
                values_target2 = get_list_of_values(exon_dic_list, target2, control_dic, nt)
                figure_creator(values_target1, values_target2, name_projects, regulation, target1,
                               target2, nt, ctrl, final_output)


if __name__ == "__main__":
    main()
