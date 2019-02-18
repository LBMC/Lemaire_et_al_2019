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
import sys
import numpy as np
# import math
# from ncephes import cprob
from rpy2.robjects.packages import importr
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import pandas as pd
import config
import stretch_calculator
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("stretch_calculator", "make_control_files_bp_ppt/"))
import exon_class
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("stretch_calculator", ""))
import group_factor
import union_dataset_function


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


# def frequency_test(obs1, tot1, obs2, tot2):
#     """
#     Proportion test.
#
#     :param obs1: (int) the count number of an amino acid X in the set of protein 1.
#     :param tot1: (int) the total number of amino acids in the set of protein 1.
#     :param obs2: (int) the count number of an amino acid X in the set of protein 2.
#     :param tot2: (int) the total number of amino acids in the set of protein 2.
#     :return: proportion test p-value
#     """
#     if obs1 == obs2:
#         return "NA"
#     if obs1 == tot1 and obs2 == tot2:
#         return "NA"
#     mean1 = float(obs1) / tot1
#     mean2 = float(obs2) / tot2
#
#     var1 = float(obs1) * (1 - mean1) * (1 - mean1) + (tot1 - obs1) * mean1 * mean1
#     var2 = float(obs2) * (1 - mean2) * (1 - mean2) + (tot2 - obs2) * mean2 * mean2
#
#     df = tot1 + tot2 - 2
#     svar = (var1 + var2) / df
#     t = (mean1 - mean2) / math.sqrt(svar * (1.0 / tot1 + 1.0 / tot2))
#     return cprob.incbet(0.5 * df, 0.5, df/(df + t * t))


def frequency_test(obs1, tot1, obs2, tot2):
    """
    Chiq test

    :param obs1: (int) the count number of an amino acid X in the set of protein 1.
    :param tot1: (int) the total number of amino acids in the set of protein 1.
    :param obs2: (int) the count number of an amino acid X in the set of protein 2.
    :param tot2: (int) the total number of amino acids in the set of protein 2.
    :return: proportion test p-value
    """

    chisq = robj.r("""

        function(vect){
            m<-matrix(vect, byrow=T, nrow=2)
            return(chisq.test(m)$p.value)
        }

                       """)
    print(obs1, tot1, obs2, tot2)
    rm1 = tot1 - obs1
    rm2 = tot2 - obs2
    vect = v.FloatVector([obs1, rm1, obs2, rm2])
    pval = float(chisq(vect)[0])
    if np.isnan(pval):
        return "NA"
    else:
        return pval


def adjust_pvalues(pvalues):
    """
    correct a list of pvalues
    :param pvalues: (list of float) list of pvalues
    :return: (list of float) list of pvalues corrected
    """
    rstats = robj.packages.importr('stats')
    pcor = np.array(rstats.p_adjust(v.FloatVector(pvalues), method="BH"))
    return list(pcor)


def write_proportion_pvalues(list_values, list_name, output, regulation, name_fig, type_fig="exon"):
    """
    Write a text file containing the pvalue of a frequency test for test sets vs control. \

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param name_fig: (string) the name of the figure
    :param type_fig: (string) the type of graphics created
    """
    list_col = ["stretch", "factor1", "factor2", "nb_bp1", "tot1",
            "nb_bp2",  "tot2", "pval"]
    data = {my_name: [] for my_name in list_col}
    filename = "%s%s_%s_exons_lvl_%s_barplot_pvalue.txt" % (output, name_fig, regulation, type_fig)
    for i in range(len(list_values) - 1):
        tot1 = len(list_values[i])
        for j in range(i + 1, len(list_values)):
            tot2 = len(list_values[j])
            for val in config.abscissa:
                nb1 = eval("sum(np.array(list_values[i]) %s)" % val)
                nb2 = eval("sum(np.array(list_values[j]) %s)" % val)
                pval = frequency_test(nb1, tot1, nb2, tot2)
                list_info = [val.replace("==", "="), list_name[i], list_name[j], nb1, tot1, nb2, tot2, pval]
                for k in range(len(list_col)):
                    data[list_col[k]].append(list_info[k])
    df = pd.DataFrame(data)
    pval_tmp = list(df["pval"])
    pval = []
    for val in pval_tmp:
        if val == "NA":
            pval.append(float("nan"))
        else:
            pval.append(val)
    pcor = adjust_pvalues(pval)
    df["pcor"] = pcor
    df = df[list_col + ["pcor"]]
    df.to_csv(filename, sep="\t", index=False)


def create_barplot(list_values, list_name, output, regulation, name_fig, type_fig="exon"):
    """
    Create a barplot showing the values of ``name_fig`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param name_fig: (string) the name of the figure
    :param type_fig: (string) the type of graphics created
    """
    list_name = [name.replace("_exons", "") for name in list_name]
    color_dic = group_factor.color_dic
    if type_fig == "factors":
        list_color = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(list_values))]
    data = []
    abscissa = [val.replace("==", "=") for val in config.abscissa]
    for i in range(len(list_values)):
        ordinate = []
        for val in config.abscissa:
            eval("ordinate.append(sum(np.array(list_values[i]) %s) / len(list_values[i]))" % val)
        data.append(go.Bar(x=abscissa, y=ordinate, name=list_name[i], marker=dict(color=color_dic[list_name[i]])))

    title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title="proportion",
            gridcolor='white',
            gridwidth=1,
            zerolinecolor='white',
            zerolinewidth=2,
        ),
        xaxis=dict(title=name_fig),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl_%s_barplot.html" % (output, name_fig, regulation, type_fig),
                    auto_open=False, validate=False)


def get_stretch_score_list(exon_list, stretch_data):
    """
    Calculate or retrieve bp and ppt score
    :param exon_list: (list of ExonClass object) List of exons
    :param stretch_data: (list of 2 int) the size and the content of the stretch wanted
    :return: (dictionary of list of int) each key of the dictionary corresponds to a particular nucleotide and \
    each list link to a nucleotide correspond to the number of stretch in each sequences  \
    before each exons in ``exon_list``.
    """
    print("Calculating the stretches of the wanted exon list")
    stretch_dic = stretch_calculator.stretch_counter(exon_list, stretch_data, config.sequence_boundaries)
    return stretch_dic


def initiate_list_of_factor(file_dir, exon_type, type_analysis):
    """
    Get the name of the exon set studied and the file containing those exon if necessary.

    :param file_dir: (string) the path where the file containing gc and AT exon set is located
    :param exon_type: (string) the type of control exon wanted
    :param type_analysis: (string) the type of analysis wanted
    :return: (2 lists):
        - list_file: (list of string or NoneType) list of the filename containing exon sets if ``type analysis = exons``
        - name_list: (list of strings) list of the name of exon sets
    """
    if type_analysis == "exon":
        at_pure_file = "%sAT_rich_exons" % file_dir
        gc_pure_file = "%sGC_rich_exons" % file_dir
        name_list = ["GC_pure_exons", "AT_pure_exons", exon_type]
        list_file = [gc_pure_file, at_pure_file, None]
    else:
        name_list = ["U1-factors", "U2-factors", exon_type]
        list_file = None
    return name_list, list_file


def extract_data(cnx, cnx_sed, list_files, list_names, pos, regulation):
    """

    :param cnx: (sqlite3 connect object) connection to fasterDB lite
    :param cnx_sed: (sqlite3 connect object) connection to sed
    :param list_files: (list of string) list of files containing exon set
    :param list_names: (list of string) the name of exon set
    :param pos: (int) the position of interest within the list ``list_files`` and ``list_names``. \
    Those 2 lists must have the same length
    :param regulation: (string) up or down
    :return: (list of ExonClass object) list of exon.
    """
    if list_files:
        exon_list = extract_exon_files(cnx, list_files[pos])
    else:
        dic_name = {"U1-factors": ["SNRPC", "SNRNP70", "DDX5_DDX17"], "U2-factors": ["U2AF2", "SF1", "SF3A3", "SF3B4"]}
        exon_list_tmp = []
        for sf_name in dic_name[list_names[pos]]:
            exon_list_tmp += union_dataset_function.get_every_events_4_a_sl(cnx_sed, sf_name, regulation)
        exon_list_tmp = union_dataset_function.washing_events_all(exon_list_tmp)
        exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list_tmp]
    print("%s : %s %s exons" % (list_names[pos], len(exon_list), regulation))
    return exon_list


def main():
    exon_class.set_debug(0)
    exon_type = "CCE"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/stretch_calculator",
                                                                      "result/stretch_calculator/")
    file_dir = os.path.realpath(os.path.dirname(__file__)).replace("src/stretch_calculator", "result/")
    if not os.path.isdir(output):
        os.mkdir(output)
    ctrl_dir = os.path.realpath(os.path.dirname(__file__)) + "/control_dictionaries/"
    sys.path.insert(0, ctrl_dir)
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src/stretch_calculator", "data/fasterDB_lite.db")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/stretch_calculator", "data/sed.db")
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    type_factors = ["exon", "spliceosome", "spliceosome"]
    regulations = ["down", "down", "up"]
    for i in range(len(type_factors)):
        type_analysis = type_factors[i]
        regulation = regulations[i]
        name_file, list_file = initiate_list_of_factor(file_dir, exon_type, type_analysis)
        dict_stretch_3ss = {"X".join(map(str, stretch_data)): {nt: [] for nt in config.nt_list} for stretch_data in config.stretches}
        for i in range(len(name_file)):
            if name_file[i] != exon_type:
                exon_list = extract_data(cnx, cnx_sed, list_file, name_file, i, regulation)
                for stretch_data in config.stretches:
                    stretch_dic = get_stretch_score_list(exon_list, stretch_data)
                    for nt in config.nt_list:
                        dict_stretch_3ss["X".join(map(str, stretch_data))][nt].append(stretch_dic[nt])
            else:
                for stretch_data in config.stretches:
                    mod = __import__("%s_stretches" % exon_type)
                    st_name = "X".join(map(str, stretch_data))
                    ctrl_dic = eval("mod.stretch_%s" % st_name)
                    for nt in config.nt_list:
                        dict_stretch_3ss[st_name][nt].append(ctrl_dic[nt])

        for stretch_data in config.stretches:
            st_name = "X".join(map(str, stretch_data))
            for nt in config.nt_list:
                create_barplot(dict_stretch_3ss[st_name][nt], name_file, output, regulation, "nb_stretch_%s-%s_%s_nt" % (stretch_data[1], stretch_data[0], nt), type_analysis)
                write_proportion_pvalues(dict_stretch_3ss[st_name][nt], name_file, output, regulation, "nb_stretch_%s-%s_%s_nt" % (stretch_data[1], stretch_data[0], nt), type_analysis)


if __name__ == "__main__":
    main()