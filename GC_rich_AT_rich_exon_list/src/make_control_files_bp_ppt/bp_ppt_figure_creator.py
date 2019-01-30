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
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from rpy2.robjects.packages import importr
import warnings
from rpy2.rinterface import RRuntimeWarning
import exon_class
import function
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("make_control_files_bp_ppt", ""))
import group_factor
import union_dataset_function
import pandas as pd

abscissa = ["<=2", ">=3"]


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
        cs1 = make_col_scheme(chars=c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M'), groups=c('g1','g2','g3','g4','g5', 'g6', 'g7', 'g8', 'g9', 'g10'),cols=c('limegreen','brown1','gold','dodgerblue3','darkorange', "brown1", "limegreen", "dodgerblue3", "darkorchid3", "dodgerblue3"), name='custom1')

        p1 = ggseqlogo(mys_seq,  method = "bit", col_scheme=cs1, namespace = c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M')) + theme_logo() + scale_x_discrete(limits = as.character(seq(1,size, by=1)), labels = as.character(seq(1,size, by=2)), breaks = as.character(seq(1, size, by=2))) + theme(axis.title.y=element_text(size=s1+25), legend.position="none")
        p1 = p1 + ggtitle(mytitle) +  theme(plot.title = element_text(hjust = 0.5))


        p1 = p1 + theme(axis.text=element_text(size=s1 + 25), plot.title = element_text(size=s1 + 30))
        p1 = p1 + scale_y_discrete(limits = c(0, 0.5, 1), labels = as.character(seq(0,1, length=3)), breaks = as.character(seq(0,1, length=3)), expand = c(0,0.05))
        #p1 = p1 + ylim(0,1)
        png(file=paste(name_file,"_weblogo.png"),height=149 * 2,width=52 * size * 2 )
        print(p1)
        dev.off()
    }
    """)
    weblogo_maker(v.StrVector(sequence_list), v.StrVector([output + sequence_name]), v.StrVector([sequence_name]), v.IntVector([len(sequence_list[0])]))


def create_figure(list_values, list_name, output, regulation, name_fig, type_fig="exon"):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of figure group or U1
    """
    filename = "%s%s_%s_exons_lvl_%s.html" % (output, name_fig, regulation, type_fig)
    list_name = [name.replace("_exons", "") for name in list_name]
    list_values[-1] = list(map(float, list_values[-1]))
    color_dic = group_factor.color_dic
    color_b_dic = group_factor.color_dic_bright
    data = []
    title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)
    with open(filename.replace("html", "txt"), "w") as outfile:
        for i in range(len(list_values) - 1):
            for j in range(i + 1, len(list_values)):
                pval = mann_withney_test_r(list_values[i], list_values[j])
                outfile.write("%s vs %s, p=%.2E\n" % (list_name[i], list_name[j], pval))

    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "fillcolor": color_b_dic[list_name[i]], "opacity": 1,
                     "line": {"color": "black"},
                     "box": {"visible": True, "fillcolor": color_dic[list_name[i]]}, "meanline": {"visible": False}})

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title=name_fig,
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
        shapes=[dict(type="line", x0=-0.5, y0=np.median(list_values[-1]), x1=len(list_values) - 0.5,
                     y1=np.median(list_values[-1]),
                     line=dict(color=color_dic[list_name[-1]]))]
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename=filename,
                        auto_open=False, validate=False)


def create_barplot_with_error_bar(list_values, list_name, output, regulation, name_fig, type_fig="exon"):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of figure group or U1
    """
    list_name = [name.replace("_exons", "") for name in list_name]
    list_values[-1] = list(map(float, list_values[-1]))
    color_dic = group_factor.color_dic
    title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)
    for i in range(len(list_values) - 1):
        pval = mann_withney_test_r(list_values[i], list_values[-1])
        title += "<br> %s vs %s, p=%.2E" % (list_name[i], list_name[-1], pval)
    if len(list_values) == 3:
        pval = mann_withney_test_r(list_values[0], list_values[1])
        title += "<br> %s vs %s, p=%.2E" % (list_name[0], list_name[1], pval)

    mean_values = [np.nanmean(values) for values in list_values]
    sd_values = [np.nanstd(values) for values in list_values]
    color_list = [color_dic[my_name] for my_name in list_name]
    trace = go.Bar(x=list_name,
                   y=mean_values,
                   marker=dict(color=color_list),
                   error_y = dict(type="data",
                                  array=sd_values,
                                  visible=True))
    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title=name_fig,
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
    plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl_%s_error_bar.html" % (output, name_fig, regulation, type_fig),
                    auto_open=False, validate=False)


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
    list_col = ["nb_branch_point", "factor1", "factor2", "nb_bp1", "tot1",
            "nb_bp2",  "tot2", "pval"]
    data = {my_name: [] for my_name in list_col}
    filename = "%s%s_%s_exons_lvl_%s_barplot_pvalue.txt" % (output, name_fig, regulation, type_fig)
    for i in range(len(list_values) - 1):
        tot1 = len(list_values[i])
        for j in range(i + 1, len(list_values)):
            tot2 = len(list_values[j])
            for val in abscissa:
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
    data = []
    abscissa2 = [val.replace("==", "=") for val in abscissa]
    for i in range(len(list_values)):
        ordinate = []
        for val in abscissa2:
            eval("ordinate.append(sum(np.array(list_values[i]) %s) / len(list_values[i]))" % val)
        data.append(go.Bar(x=abscissa2, y=ordinate, name=list_name[i], marker=dict(color=color_dic[list_name[i]])))

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


def mann_withney_test_r(list_values1, list_values2):
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def get_bp_ppt_score_list(output, exon_list, name_list, size, regulation):
    """
    Calculate or retrieve bp and ppt score
    :param output: (string) path where the result will be stored or retrieved
    :param exon_list: (list of ExonClass object) List of exons
    :param name_list: (string) the name of the exons list
    :param size: (int) the size of window wanted
    :param regulation: (string) up or down
    :return: (2 list of floats) list of bp score and list of ppt score
    """
    name_store_file = "%s%s_%s_%s_bp_ppt_score.py" % (output, name_list, size, regulation)
    print(name_store_file)
    if not os.path.isfile(name_store_file):
        print("Calculating ppt bp score using SVM BP FINDER")
        bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, bp_seq_list, ag_count_list, hbound_list, uaa_list, una_list = function.bp_ppt_calculator(exon_list, size)
        with open(name_store_file, "w") as outfile:
            outfile.write("bp_score_list=%s\n" % str(bp_score_list))
            outfile.write("ppt_score_list=%s\n" % str(ppt_score_list))
            outfile.write("nb_bp_list=%s\n" % str(nb_bp_list))
            outfile.write("nb_good_bp_list=%s\n" % str(nb_good_bp_list))
            outfile.write("bp_seq=%s\n" % str(bp_seq_list))
            outfile.write("ag_count=%s\n" % str(ag_count_list))
            outfile.write("hbound=%s\n" % str(hbound_list))
            outfile.write("uaa_count=%s\n" % str(uaa_list))
            outfile.write("una_count=%s\n" % str(una_list))
    else:
        print("recovering ppt and bp score already stored in %s" % name_store_file)
        sys.path.insert(0, output)
        stored = __import__("%s_%s_%s_bp_ppt_score" % (name_list, size, regulation))
        bp_score_list = stored.bp_score_list
        ppt_score_list = stored.ppt_score_list
        nb_bp_list = stored.nb_bp_list
        nb_good_bp_list = stored.nb_good_bp_list
        bp_seq_list = stored.bp_seq
        ag_count_list = stored.ag_count
        hbound_list = stored.hbound
        uaa_list = stored.uaa_count
        una_list = stored.una_count
    return bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, bp_seq_list, ag_count_list, hbound_list, uaa_list, una_list


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
        name_list = ["SNRNP70", "SNRPC", "U2AF2", "SF1", exon_type]
        list_file = None
    return name_list, list_file


def extract_data(cnx, cnx_sed, list_files, list_names, pos, regulation="down"):
    """

    :param cnx: (sqlite3 connect object) connection to fasterDB lite
    :param cnx_sed: (sqlite3 connect object) connection to sed
    :param list_files: (list of string) list of files containing exon set
    :param list_names: (list of string) the name of exon set
    :param pos: (int) the position of interest within the list ``list_files`` and ``list_names``. \
    Those 2 lists must hace the same lenght
    :param regulation: (string) up or down
    :return: (list of ExonClass object) list of exon.
    """
    if list_files:
        exon_list = extract_exon_files(cnx, list_files[pos])
    else:
        exon_list_tmp = union_dataset_function.get_every_events_4_a_sl(cnx_sed, list_names[pos], regulation)
        exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list_tmp]
    print("%s : %s %s exons" % (list_names[pos], len(exon_list), regulation))
    return exon_list

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
    type_factors = ["exon", "spliceosome", "spliceosome"]
    regulations = ["down", "down", "up"]

    for i in range(len(type_factors)):
        type_analysis = type_factors[i]
        regulation = regulations[i]
        name_file, list_file = initiate_list_of_factor(file_dir, exon_type, type_analysis)
        dict_score_3ss = {i:{"bp_score_list":[], "ppt_score_list":[], "nb_bp_list":[], "nb_good_bp_list":[], "ag_count": [], "hbound": [], "uaa_count": [], "una_count": []} for i in [100, 50, 35, 25]}
        list_force_acceptor = []
        list_force_donor = []
        for i in range(len(name_file)):
            if name_file[i] != exon_type:
                exon_list = extract_data(cnx, cnx_sed, list_file, name_file, i, regulation)
                list_force_acceptor.append(get_redundant_list_of_value(cnx_sed, exon_list, "force_acceptor"))
                list_force_donor.append(get_redundant_list_of_value(cnx_sed, exon_list, "force_donor"))
                for size in dict_score_3ss.keys():
                    bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, bp_seq_list, ag_count_list, hbound_list, uaa_list, una_list = get_bp_ppt_score_list(ctrl_output, exon_list, name_file[i], size, regulation)
                    dict_score_3ss[size]["bp_score_list"].append(bp_score_list)
                    dict_score_3ss[size]["ppt_score_list"].append(ppt_score_list)
                    dict_score_3ss[size]["nb_bp_list"].append(nb_bp_list)
                    dict_score_3ss[size]["nb_good_bp_list"].append(nb_good_bp_list)
                    web_logo_creator(bp_seq_list, "%s_%s_exons_%s_nt" % (name_file[i], regulation, size), output)
                    dict_score_3ss[size]["ag_count"].append(ag_count_list)
                    dict_score_3ss[size]["hbound"].append(hbound_list)
                    dict_score_3ss[size]["uaa_count"].append(uaa_list)
                    dict_score_3ss[size]["una_count"].append(una_list)


            else:
                for size in dict_score_3ss.keys():
                    mod = __import__("%s_%s_bp_ppt_score" % (exon_type, size))
                    dict_score_3ss[size]["bp_score_list"].append(mod.bp_score)
                    dict_score_3ss[size]["ppt_score_list"].append(mod.ppt_score)
                    dict_score_3ss[size]["nb_bp_list"].append(mod.nb_bp)
                    dict_score_3ss[size]["nb_good_bp_list"].append(mod.nb_good_bp)
                    web_logo_creator(mod.bp_seq, "%s_%s_exons_%s_nt" % (name_file[i], regulation, size), output)
                    dict_score_3ss[size]["ag_count"].append(mod.ag_count)
                    dict_score_3ss[size]["hbound"].append(mod.hbound)
                    dict_score_3ss[size]["uaa_count"].append(mod.uaa_count)
                    dict_score_3ss[size]["una_count"].append(mod.una_count)
                list_force_acceptor.append(get_control_exon_information(cnx_sed, exon_type, "force_acceptor"))
                list_force_donor.append(get_control_exon_information(cnx_sed, exon_type, "force_donor"))

        for size in dict_score_3ss.keys():
            print("------------> %s nt " % size)
            create_barplot(dict_score_3ss[size]["nb_bp_list"], name_file, output, regulation, "prop_nb_branch_point_(%snt)" % size, type_analysis)
            create_barplot(dict_score_3ss[size]["nb_good_bp_list"], name_file, output, regulation, "prop_nb_good_branch_point_(%snt)" % size, type_analysis)
            write_proportion_pvalues(dict_score_3ss[size]["nb_good_bp_list"], name_file, output, regulation, "prop_nb_good_branch_point_(%snt)" % size, type_analysis)
            create_barplot(dict_score_3ss[size]["ag_count"], name_file, output, regulation,
                           "AG_count_downstream_bp(%snt)" % size, type_analysis)
            write_proportion_pvalues(dict_score_3ss[size]["ag_count"], name_file, output, regulation,
                                     "AG_count_downstream_bp(%snt)" % size, type_analysis)

            create_barplot(dict_score_3ss[size]["uaa_count"], name_file, output, regulation,
                           "UAA_count(%snt)" % size, type_analysis)
            write_proportion_pvalues(dict_score_3ss[size]["uaa_count"], name_file, output, regulation,
                                     "UAA_count(%snt)" % size, type_analysis)
            create_barplot(dict_score_3ss[size]["una_count"], name_file, output, regulation,
                           "UNA_count(%snt)" % size, type_analysis)
            write_proportion_pvalues(dict_score_3ss[size]["una_count"], name_file, output, regulation,
                                     "UNA_count(%snt)" % size, type_analysis)
            create_figure(dict_score_3ss[size]["hbound"], name_file, output, regulation, "nb_h_bound_%s_nt" % size, type_analysis)
            create_barplot_with_error_bar(dict_score_3ss[size]["hbound"], name_file, output, regulation, "nb_h_bound_%s_nt" % size,
                          type_analysis)
        create_figure(list_force_acceptor, name_file, output, regulation, "force_acceptor", type_analysis)
        create_barplot_with_error_bar(list_force_acceptor, name_file, output, regulation, "force_acceptor", type_analysis)
        create_figure(list_force_donor, name_file, output, regulation, "force_donor", type_analysis)
        create_barplot_with_error_bar(list_force_donor, name_file, output, regulation, "force_donor", type_analysis)



if __name__ == "__main__":
    main()