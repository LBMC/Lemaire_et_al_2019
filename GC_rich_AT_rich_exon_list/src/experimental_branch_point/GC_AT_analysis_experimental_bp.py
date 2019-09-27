#!/usr/bin/env python3

"""
Description:
    The goal of this script is to count the number of experimental \
    branch point in the intron upstream GC and AT exons
"""

import os
import sqlite3
import pandas as pd
import sys
import plotly.graph_objs as go
import plotly
import numpy as np
import lazyparser as lp
base = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, base)
import union_dataset_function as udf
import group_factor
sys.path.insert(0, base + "/make_control_files_bp_ppt")
import bp_ppt_figure_creator as bpfc
abscissa = ["<=2", ">2"]


def get_ctrl_exons(cnx, exon_type, exon2remove):
    """
    Get CCE exons

    :param cnx: (sqlite3 object) allow connection to fasterdb database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of list of 2int) list of exon to remove frome the control list of exons
    :return: (list of list of 2 int) list of CCE exons
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT id_gene, pos_on_gene
                   FROM exons
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT id_gene, pos_on_gene
                   FROM exons
                   AND id_gene = t2.id
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    print("number of control exon before removing bad ones : %s" % len(result))
    nresult = [list(exon) for exon in result if [exon[0], exon[1]] not in exon2remove]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    return nresult


def read_file(filename):
    """
    Read a file full of exons.

    :param filename: (str) a file containing exons
    :return: (list of str) list of exons
    """
    list_exons = []
    with open(filename, "r") as infile:
        for line in infile:
            line = line.replace("\n", "")
            line = line.split("\t")
            line = [int(val) if val.isdigit() else val for val in line]
            list_exons.append(line)
    return list_exons


def get_intron_coordinates(cnx, list_exons, type_exons):
    """

    :param cnx: (sqlite3 connectio object) connection to fasterdb lite
    :param list_exons: (list of 2 int) list of GC and AT exons
    :param type_exons: (list of str) type of the exon GC or AT
    :return:
    """
    dic = {1: "+", -1: "-"}
    count = 0
    cursor = cnx.cursor()
    query = """
    SELECT t1.chromosome, t2.start_on_chromosome - 1, t2.end_on_chromosome, 
           t1.strand
    FROM genes t1, introns t2
    WHERE t1.id = t2.id_gene
    AND   t2.id_gene = %s
    AND   t2.pos_on_gene = %s
    """
    content_intron_bed = []
    for exon, type_exon in zip(list_exons, type_exons):
        cursor.execute(query % (exon[0], exon[1] - 1))
        result = cursor.fetchall()
        if len(result) > 1:
            raise IndexError("multiple introns found : Intron %s" % exon)
        if len(result) == 1:
            intron = result[0]
            intron = list(intron)
            if intron[1] < intron[2]:
                intron[3] = dic[intron[3]]
                if intron[3] == "+" and intron[2] - intron[1] > 100:
                    intron[1] = intron[2] - 100
                if intron[3] == "-" and intron[2] - intron[1] > 100:
                    intron[2] = intron[1] + 100
                content_intron_bed.append([intron[0], intron[1], intron[2]] +
                                          ["%s_%s" % (exon[0],exon[1] - 1)] +
                                          [".", intron[3], type_exon])
            else:
                count += 1
    print("Exon with a negative or nul size : %s" % count)
    return content_intron_bed


def intersection_exon(intron, list_bp):
    """
    Get the number of branch points entirely within the current intron.

    :param intron: (list of data) list of data corresponding to a TAD \
    in a bed order.
    :param list_bp: (list of list of data) a list of line in a bed file \
    containing branch points
    :return: (list of str) list of exon within the current tad
    """
    intron_start = intron[1]
    intron_stop = intron[2]
    if intron_stop <= intron_start:
        raise ValueError("The size of the intron %s if negative or nul" %
                         intron[3])
    res = []
    for line in list_bp:
        if line[2] > line[1]:
            if  line[0].replace("chr", "") == intron[0] and \
                    line[5] == intron[5] and \
                    (intron_start <= line[1] <= intron_stop and
                     intron_start <= line[2] <= intron_stop):
                res.append(line)
    return len(res)


def get_intron_bp_data(list_intron, list_bp):
    """
    Get the number of bp within each intron.

    :param list_intron: (list of list of values) list of introns
    :param list_bp: (list of list of value) the data stored in a bed file
    :return: (pandas dataframe)
    """
    print("getting intron data")
    dic = {"intron": [], "type": [], "nb_bp": []}
    tot = len(list_intron)
    count = 0
    for intron in list_intron:
        nb_bp = intersection_exon(intron, list_bp)
        dic["intron"].append(intron[3])
        dic["type"].append(intron[6])
        dic["nb_bp"].append(nb_bp)
        count += 1
        sys.stdout.write("%s/%s       \r" % (count, tot))
    df = pd.DataFrame(dic)
    return df


def create_barplot(df, output, name_fig):
    """
    Create a barplot showing the values of ``name_fig`` for every exon list regulated by ``list_name`` factor

    :param df: (pandas dataframe) a pandas dataframe containing the number of experimental branchpoint
    :param regulation: (string) up or down
    :param name_fig: (string) the name of the figure
    """
    list_name = list(np.unique(df["type"]))
    list_name = [list_name[-1]] + list_name[:-1]
    print(list_name)
    list_values = [df[(df["type"] == mtype)]["nb_bp"].values
                   for mtype in list_name]
    color_dic = group_factor.color_dic
    default_colors = [color_dic["GC_pure"], color_dic["AT_pure"], "#F3A431",
                      "#A000A0", "#5c0001", color_dic["CCE"]]
    if len(default_colors) != len(list_values):
        default_colors = default_colors[0:len(list_values) - 1] + \
                         [default_colors[-1]]
    data = []
    abscissa2 = [val.replace("==", "=") for val in abscissa]
    for i in range(len(list_values)):
        try:
            if "GC" in list_name[i] or "AT" in list_name[i]:
                mcol = color_dic[list_name[i].replace("-exons", "_pure")]
            else:
                mcol = color_dic[list_name[i].replace("-exons", "")]
        except KeyError:
            mcol = default_colors[i]
        ordinate = []
        for val in abscissa:
            eval("ordinate.append(sum(list_values[i] %s) / len(list_values[i]))" % val)
        data.append(go.Bar(x=abscissa2, y=ordinate, name=list_name[i], marker=dict(color=mcol)))

    title = """%s of exons regulated by different factors""" % name_fig

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            # autotick=True,
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
    filename = "%s/%s_barplot.html" % (output, name_fig)
    plotly.offline.plot(fig, filename=filename,
                        auto_open=False, validate=False)
    write_proportion_pvalues(list_values, list_name,
                             filename.replace(".html", "_stat.txt"), name_fig)


def write_proportion_pvalues(list_values, list_name, filename, name_fig):
    """
    Write a text file containing the pvalue of a frequency test for test sets vs control. \

    :param list_values: (list of list of float) each sublist corresponds \
    to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different group \
    of intron studied.
    :param filename: (str) the name of the file to create
    :param name_fig: (string) the name of the figure
    """
    list_col = [name_fig, "factor1", "factor2", "nb_bp1", "tot1",
                "nb_bp2",  "tot2", "pval"]
    data = {my_name: [] for my_name in list_col}
    for i in range(len(list_values) - 1):
        tot1 = len(list_values[i])
        for j in range(i + 1, len(list_values)):
            tot2 = len(list_values[j])
            for val in abscissa:
                nb1 = eval("sum(list_values[i] %s)" % val)
                nb2 = eval("sum(list_values[j] %s)" % val)
                pval = bpfc.frequency_test(nb1, tot1, nb2, tot2)
                list_info = [val.replace("==", "="), list_name[i],
                             list_name[j], nb1, tot1, nb2, tot2, pval]
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
    pcor = bpfc.adjust_pvalues(pval)
    df["pcor"] = pcor
    df = df[list_col + ["pcor"]]
    df.to_csv(filename, sep="\t", index=False)



@lp.parse(branch_point_file="file", at_file="file", gc_file="file",
          fasterdb="file", sed="file", exon_type=["CCE"], output="dir")
def main(branch_point_file, name_bp_file, at_file, gc_file, fasterdb, seddb,
         output, exon_type="CCE"):
    """
    Create a GC/AT barplots with experimental branch points.

    :param branch_point_file: (str) a bed file containing branch point
    :param name_bp_file: (str) the name of the experimental bp file
    :param at_file: (str) a file containing AT exons
    :param gc_file: (str) a file containing GC exons
    :param fasterdb: (str) path to fasterdb database
    :param seddb: (str) path to sed database
    :param output: (str) folder where the figures will be created
    :param exon_type: (str) the type of control exons
    """

    result_file = "%s/intron_experimental_%s_table.txt" % (output,
                                                           name_bp_file)
    if not os.path.isfile(result_file):
        at_exon = read_file(at_file)
        gc_exon = read_file(gc_file)
        list_bp = read_file(branch_point_file)
        print(len(list_bp))
        cnx = sqlite3.connect(fasterdb)
        cnx_sed = sqlite3.connect(seddb)
        exon2remove = [list(map(int, exon))
                       for exon in udf.get_exon_regulated_by_sf(cnx_sed,
                                                                "down")]
        ctrl_exons = get_ctrl_exons(cnx, exon_type, exon2remove)
        exon_list = gc_exon + at_exon + ctrl_exons
        type_exon = ["GC-exons"] * len(gc_exon) + \
                    ["AT-exons"] * len(at_exon) + \
                    ["%s-exons" % exon_type] * len(ctrl_exons)
        intron_data = get_intron_coordinates(cnx, exon_list, type_exon)
        df = get_intron_bp_data(intron_data, list_bp)
        print(df.head())
        df.to_csv(result_file, sep="\t",
                  index=False)
        cnx.close()
        cnx_sed.close()
    else:
        print("Recovering %s" % result_file)
        df = pd.read_csv(result_file, sep="\t")
    create_barplot(df, output, os.path.basename(result_file).replace(".txt", ""))


if __name__ == "__main__":
    main()


