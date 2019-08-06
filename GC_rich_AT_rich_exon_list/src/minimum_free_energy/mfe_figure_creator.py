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
import math
import stat_mfe
import seaborn as sns
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from rpy2.robjects.packages import importr
import pandas as pd
import warnings
from rpy2.rinterface import RRuntimeWarning
import exon_class
import function
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("minimum_free_energy", ""))
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


def web_logo_creator(sequence_list, sequence_name, output):
    """
    :param sequence_list: (tuple of strings) - list of sequences
    :param sequence_name: (string) name of the sequence
    :param output: (string) the folder where the results will be created
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    weblogo_maker = robj.r("""
    library("ggplot2")
    library("ggseqlogo")

    function(mys_seq, name_file, mytitle, size){
        s1 = 15
        cs1 = make_col_scheme(chars=c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M'), 
                              groups=c('g1','g2','g3','g4','g5', 'g6', 'g7', 'g8', 'g9', 'g10'),
                              cols=c('limegreen','brown1','gold','dodgerblue3','darkorange', "brown1", 
                                     "limegreen", "dodgerblue3", "darkorchid3", "dodgerblue3"), 
                              name='custom1')

        p1 = ggseqlogo(mys_seq,  method = "bit", col_scheme=cs1, 
                       namespace = c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M')) + 
                       theme_logo() + 
                       scale_x_discrete(limits = as.character(seq(1,size, by=1)), 
                       labels = as.character(seq(1,size, by=2)), breaks = as.character(seq(1, size, by=2))) + 
                       theme(axis.title.y=element_text(size=s1+25), legend.position="none")
        p1 = p1 + ggtitle(mytitle) +  theme(plot.title = element_text(hjust = 0.5))


        p1 = p1 + theme(axis.text=element_text(size=s1 + 25), plot.title = element_text(size=s1 + 30))
        p1 = p1 + scale_y_discrete(limits = c(0, 0.5, 1), labels = as.character(seq(0,1, length=3)), 
                                   breaks = as.character(seq(0,1, length=3)), expand = c(0,0.05))
        #p1 = p1 + ylim(0,1)
        png(file=paste(name_file,"_weblogo.png"),height=149 * 2,width=52 * size * 2 )
        print(p1)
        dev.off()
    }
    """)
    weblogo_maker(v.StrVector(sequence_list), v.StrVector([output + sequence_name]),
                  v.StrVector([sequence_name]), v.IntVector([len(sequence_list[0])]))


def create_figure(list_values, list_name, output, regulation, splicing_site, type_fig,
                  num_fig=""):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param splicing_site: (string) 5ss or 3SS
    :param type_fig: (string) exon or spliceosome
    """
    list_name = [name.replace("_exons", "") for name in list_name]
    color_dic = group_factor.color_dic
    color_b_dic = group_factor.color_dic_bright
    data = []
    title = """Minimum free energy of %s of %s exon regulated by different factors""" % (splicing_site, regulation)
    default_colors = [color_dic["GC_pure"], color_dic["AT_pure"], color_dic["CCE"]]
    default_bright = [color_b_dic["GC_pure"], color_b_dic["AT_pure"], color_b_dic["CCE"]]
    for i in range(len(list_values)):
        try:
            cur_color = color_dic[list_name[i]]
            color_b = color_b_dic[list_name[i]]
        except KeyError:
            cur_color = default_colors[i]
            color_b = default_bright[i]
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "fillcolor": color_b, "opacity": 1,
                     "line": {"color": "black"},
                     "box": {"visible": True, "fillcolor": cur_color}, "meanline": {"visible": False}})

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title="MFE",
            gridcolor='white',
            gridwidth=1,
            zerolinecolor='white',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        showlegend=True,
        # shapes=[dict(type="line", x0=-0.5, y0=np.median(list_values[-1]), x1=len(list_values) - 0.5,
        #              y1=np.median(list_values[-1]),
        #              line=dict(color=color_dic[list_name[-1]]))]
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%sMFE_%s_%s.html" % (output, num_fig, splicing_site, type_fig),
                        auto_open=False, validate=False)


def create_figure_error_bar(list_values, list_name, output, regulation, splicing_site, type_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param splicing_site: (string) 5ss or 3SS
    :param type_fig: (string) exon or spliceosome
    """
    list_name = [name.replace("_exons", "") for name in list_name]
    color_dic = group_factor.color_dic
    title = """Minimum free energy of %s of %s exon regulated by different factors""" % (splicing_site, regulation)

    mean_values = [np.nanmean(values) for values in list_values]
    sd_values = [np.nanstd(values) for values in list_values]
    color_list = [color_dic[my_name] for my_name in list_name]
    trace = go.Bar(x=list_name,
                   y=mean_values,
                   marker=dict(color=color_list),
                   error_y=dict(type="data",
                                array=sd_values,
                                visible=True))

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title="MFE",
            gridcolor='white',
            gridwidth=1,
            zerolinecolor='white',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='white',
        plot_bgcolor='white'
    )

    fig = {"data": [trace], "layout": layout}
    plotly.offline.plot(fig, filename="%sMFE_%s_%s_error_bar.html" % (output, splicing_site, type_fig),
                        auto_open=False, validate=False)


def adjust_pvalues(pvalues):
    """
    correct a list of pvalues
    :param pvalues: (list of float) list of pvalues
    :return: (list of float) list of pvalues corrected
    """
    rstats = robj.packages.importr('stats')
    pcor = np.array(rstats.p_adjust(v.FloatVector(pvalues), method="BH"))
    return list(pcor)


def write_proportion_pvalues(list_values, list_name, output, splicing_site, type_fig="exon"):
    """
    Write a text file containing the pvalue of a frequency test for test sets vs control. \

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param splicing_site: (string) 5ss or 3ss
    :param type_fig: (string) the type of graphics created
    """
    data = {"factor1": [], "factor2": [], "mean1": [], "mean2": [], "pval": []}
    filename = "%sMFE_%s_%s.txt" % (output, splicing_site, type_fig)
    for i in range(len(list_values) - 1):
        for j in range(i + 1, len(list_values)):
            pval = mann_withney_test_r(list_values[i], list_values[j])
            data["factor1"].append(list_name[i])
            data["factor2"].append(list_name[j])
            data["mean1"].append(np.mean(list_values[i]))
            data["mean2"].append(np.mean(list_values[j]))
            data["pval"].append(pval)
    pcor = adjust_pvalues(data["pval"])
    data["pcor"] = pcor
    df = pd.DataFrame(data)
    df = df[["factor1", "factor2", "mean1", "mean2", "pval", "pcor"]]
    df.to_csv(filename, sep="\t", index=False)


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def get_mfe_score_list(output, exon_list, name_list):
    """
    Calculate or retrieve mfe score
    :param output: (string) path where the result will be stored or retrieved
    :param exon_list: (list of ExonClass object) List of exons
    :param name_list: (string) the name of the exons list
    :return: (2 list of floats) list of bp score and list of ppt score
    """
    name_store_file = "%s%s_mfe_score.py" % (output, name_list)
    print(name_store_file)
    if not os.path.isfile(name_store_file):
        print("Calculating mfe score using RNAfold")
        mfe_3ss, mfe_5ss = function.mfe_calculator(exon_list)
        with open(name_store_file, "w") as outfile:
            outfile.write("mfe_3ss=" + str(mfe_3ss) + "\n")
            outfile.write("mfe_5ss=" + str(mfe_5ss) + "\n")
    else:
        print("recovering mfe score already stored in %s" % name_store_file)
        sys.path.insert(0, output)
        stored = __import__("%s_mfe_score" % name_list)
        mfe_3ss = stored.mfe_3ss
        mfe_5ss = stored.mfe_5ss
    return mfe_3ss, mfe_5ss


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
        name_list = ["DDX5_DDX17", "SNRNP70", "SNRPC", "U2AF1", "U2AF2", "SF1", "SF3A3", "SF3B4", exon_type]
        list_file = None
    return name_list, list_file


def extract_data(cnx, cnx_sed, list_files, list_names, pos):
    """

    :param cnx: (sqlite3 connect object) connection to fasterDB lite
    :param cnx_sed: (sqlite3 connect object) connection to sed
    :param list_files: (list of string) list of files containing exon set
    :param list_names: (list of string) the name of exon set
    :param pos: (int) the position of interest within the list ``list_files`` and ``list_names``. \
    Those 2 lists must hace the same lenght
    :return: (list of ExonClass object) list of exon.
    """
    if list_files:
        exon_list = extract_exon_files(cnx, list_files[pos])
    else:
        exon_list_tmp = union_dataset_function.get_every_events_4_a_sl(cnx_sed, list_names[pos], "down")
        exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list_tmp]
    print("%s : %s exons" % (list_names[pos], len(exon_list)))
    return exon_list


def dataframe_creator(list_values, list_name, output, regulation, name_df, type_fig, fig_num=""):
    """
    Create a pandas dataframe.

    :param list_values: (list of list of float) list of values
    :param list_name: (list of string) the names of the exon list used to get each sublist on ``list_values`` objects.
    :param name_df: (string) the type of values displayed in ``list_values``
    :param output: (string) folder where the result will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of graphic build
    :param fig_num: (str) a number
    :return: (pandas Dataframe) a dataFrame with the values in ``list_values`` and the names in ``list_names``
    """
    sns.set()
    sns.set_context("poster")
    new_values = np.hstack(list_values)
    new_names = np.hstack([[list_name[i]] * len(list_values[i]) for i in range(len(list_values))])
    new_values = new_values.astype(np.float)
    new_names = list(new_names[~np.isnan(new_values)])
    new_values = list(new_values[~np.isnan(new_values)])
    df = pd.DataFrame({"values": new_values, "project": new_names})
    # Transformation of Anscombe
    df["values"] = list(map(math.sqrt, (df["values"].values * -1) + 3/8))
    g = sns.FacetGrid(data=df, row="project", height=20)
    g.map(sns.distplot, "values")
    filename = "%s%s%s_%s_dataframe_%s_table.txt" % (output, fig_num, name_df, regulation, type_fig)
    g.savefig(filename.replace("table.txt", "displot.pdf"), format="pdf")
    df.to_csv(filename, index=False, sep="\t")
    if type_fig == "exon":
        new_df = stat_mfe.anova_nt_stats(df, filename)
    else:
        new_df = stat_mfe.anova_nt_stats_spliceosome(df, filename)
    new_df.to_csv(filename.replace("table.txt", "stat.txt"), index=False, sep="\t")


def main():
    exon_class.set_debug(0)
    exon_type = "CCE"
    ctrl_output = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy",
                                                                      "result/minimum_free_energy/")
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy",
                                                                 "result/minimum_free_energy/")
    file_dir = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy", "result/")
    if not os.path.isdir(ctrl_output):
        os.mkdir(ctrl_output)
    if not os.path.isdir(output):
        os.mkdir(output)
    ctrl_dir = os.path.realpath(os.path.dirname(__file__)) + "/control_dictionaries/"
    sys.path.insert(0, ctrl_dir)
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy", "data/fasterDB_lite.db")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy", "data/sed.db")
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    type_factors = ["exon", "spliceosome"]
    for type_analysis in type_factors:
        name_file, list_file = initiate_list_of_factor(file_dir, exon_type, type_analysis)
        mfe_3ss_score = []
        mfe_5ss_score = []
        for i in range(len(name_file)):
            if name_file[i] != exon_type:
                exon_list = extract_data(cnx, cnx_sed, list_file, name_file, i)
                mfe_3ss, mfe_5ss = get_mfe_score_list(ctrl_output, exon_list, name_file[i])
                mfe_3ss_score.append(mfe_3ss)
                mfe_5ss_score.append(mfe_5ss)
            else:
                mod = __import__("%s_mfe" % exon_type)
                mfe_3ss_score.append(mod.mfe_3ss)
                mfe_5ss_score.append(mod.mfe_5ss)

        create_figure(mfe_3ss_score, name_file, output, "down", "3SS", type_analysis)
        dataframe_creator(mfe_3ss_score, name_file, output, "down", "3SS", type_analysis)
        # create_figure_error_bar(mfe_3ss_score, name_file, output, "down", "3SS", type_analysis)
        # write_proportion_pvalues(mfe_3ss_score, name_file, output, "3SS", type_analysis)
        create_figure(mfe_5ss_score, name_file, output, "down", "5SS", type_analysis)
        dataframe_creator(mfe_5ss_score, name_file, output, "down", "5SS", type_analysis)
        # create_figure_error_bar(mfe_5ss_score, name_file, output, "down", "5SS", type_analysis)
        # write_proportion_pvalues(mfe_5ss_score, name_file, output, "5SS", type_analysis)


def main_2d(list_file, name_file, exon_type, output, seddb, fasterdb,
            fig_nums=("2.1D", "2.2D")):
    """
    Create the figure 2.1D and 2.2D of the article with a given list of exons.

    :param list_file: (list of str) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (list of str) the name of each files of exons \
    given in ``list_file``
    :param exon_type: (str) the control exons
    :param output: (str) folder where the result will be created
    :param seddb: (str) path to sed database
    :param fasterdb: (str) path to fasterdb database
    :param fig_nums: (list of str) list of figure names
    :return:
    """
    exon_class.set_debug(0)
    list_file.append(None)
    name_file.append(exon_type)

    ctrl_output = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy",
                                                                      "result/minimum_free_energy/")
    if not os.path.isdir(ctrl_output):
        os.mkdir(ctrl_output)
    ctrl_dir = os.path.realpath(os.path.dirname(__file__)) + \
               "/control_dictionaries/"
    sys.path.insert(0, ctrl_dir)
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    type_analysis = "exon"

    mfe_3ss_score = []
    mfe_5ss_score = []
    for i in range(len(name_file)):
        if name_file[i] != exon_type:
            exon_list = extract_data(cnx, cnx_sed, list_file, name_file, i)
            mfe_3ss, mfe_5ss = get_mfe_score_list(ctrl_output, exon_list, name_file[i])
            mfe_3ss_score.append(mfe_3ss)
            mfe_5ss_score.append(mfe_5ss)
        else:
            mod = __import__("%s_mfe" % exon_type)
            mfe_3ss_score.append(mod.mfe_3ss)
            mfe_5ss_score.append(mod.mfe_5ss)

    create_figure(mfe_5ss_score, name_file, output, "down", "5SS",
                  type_analysis, fig_nums[0])
    dataframe_creator(mfe_5ss_score, name_file, output, "down", "5SS",
                      type_analysis, fig_nums[0])
    create_figure(mfe_3ss_score, name_file, output, "down", "3SS",
                  type_analysis, fig_nums[1])
    dataframe_creator(mfe_3ss_score, name_file, output, "down", "3SS",
                      type_analysis, fig_nums[1])
    cnx.close()
    cnx_sed.close()

if __name__ == "__main__":
    main()
