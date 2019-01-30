#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is to create from a result file produced by the script \
    ``spliceosome_regulation_enrichment`` 4 graphics showing if every factors of the spliceosome regulates more often \
    exons from the AT group or exon from the GC group.
"""

import math
import pandas as pd
import plotly
import plotly.graph_objs as go
import os
import sys
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("/GC_AT_group_regulated_U1_U2", ""))
import group_factor


def file_reader(filename):
    """
    Load the results in ``filename`` (produced by the script  ``spliceosome_regulation_enrichment``)
    :param filename: (string) a result file created by the script ``spliceosome_regulation_enrichment``
    :return: (a dictionary of  list float + string) link each comparison to its pvalue + the direction of the comparison
    """
    dic_res = {}
    with open(filename, "r") as resfile:
        for line in resfile:
            line = line.replace("\n", "")
            if "#" not in line:
                line = line.split("\t")
                val1 = float(line[1].split("/")[0])
                name1 = line[0].split("-")[0]
                val2 = float(line[3].split("/")[0])
                name2 = line[0].split("-")[1]
                if val1 > val2:
                    dic_res[line[0]] = [float(line[7]), name1]
                elif val1 < val2:
                    dic_res[line[0]] = [float(line[7]), name2]
                else:
                    dic_res[line[0]] = [float(line[7]), "NA"]
    return dic_res


def subsample_result_dic(res_dic, group_exon, group_factor):
    """
    Do a subsample of ``res_dic`` based on ``group_exon`` and ``group_factor``/

    :param res_dic: (a dictionary of  list float + string) links each comparison to its pvalue + \
    the direction of the comparison
    :param group_exon: (string) all or pure the purity of the exon set of interest
    :param group_factor: (string) all or pure the purity of the factor of interest
    :return: (a dictionary of  list float + string) subsampled dic
    """
    sub_dic = {}
    for key in res_dic.keys():
        name = key.split("-")
        if group_exon == name[2] and group_factor == name[5]:
            sub_dic[key] = res_dic[key]
    return sub_dic


def value_adapter(sub_dic):
    """
    Transforms the values in sub dict into values that will be display in the graph
    :param sub_dic: (a dictionary of  list float + string) links each comparison to its pvalue + \
    the direction of the comparison
    :return: (list of 3 lists (of float and  2 string respectively)
    """
    list_values = []
    list_factor_name = []
    list_reg = []
    dic_reg = {"NA": 1}
    neg_reg = "..."
    for key in sub_dic.keys():
        if sub_dic[key][1] not in list_reg and sub_dic[key][1] != "NA":
            list_reg.append(sub_dic[key][1])
    print(list_reg)
    if len(list_reg) > 2:
        print("error : only 2 regulations must be in the subsampled dic!")
        exit(1)
    for val in list_reg:
        dic_reg[val] = 1

    if "GC" in list_reg:
        dic_reg["GC"] = -1
        pos_reg = "GC"
        for val in list_reg:
            if val != "GC":
                neg_reg = val
    elif "AT" in list_reg:
        dic_reg["AT"] = -1
        pos_reg = "AT"
        for val in list_reg:
            if val != "AT":
                neg_reg = val
    else:
        dic_reg[list_reg[0]] = -1
        pos_reg = list_reg[0]
        neg_reg = list_reg[1]
    for key in sub_dic.keys():
        list_values.append(math.log10(sub_dic[key][0] + 0.00001) * dic_reg[sub_dic[key][1]])
        list_factor_name.append(key.split('-')[4])
    return list_values, list_factor_name, pos_reg, neg_reg


def sort_values(value_list, name_list):
    """
    Sort the values in ``value_list``
    :param value_list: (list of float) list value
    :param name_list:  (list of string) name associated to ``value_list``
    :return: (2 lists of float and string respectively) the list oreder by ascending value in ``value_list``
    """
    df = pd.DataFrame(value_list, index=name_list)
    df = df.sort_values(0)
    value_list = df.values
    value_list = [val[0] for val in value_list]
    name_list = list(df.index)
    return value_list, name_list


def figure_maker(value_list, name_list, pos_reg, neg_reg, group_exon, group_factor_content, output):
    """
    Create a barplot from the result file produced by the script ``spliceosome_regulation_enrichment`.

    :param value_list: (list of float) list of value to display in the barplot
    :param name_list: (list of string) name of list associated to each value in ``value_list``
    :param pos_reg: (string) the positive regulation
    :param neg_reg:(string) the negative regulation
    :param group_exon: (string) the type of exon sets regulated (pure or all)
    :param group_factor_content: (string) the type of exon sets regulated by a factor(pure or all)
    :param output: (string) path where te result will be created
    """
    for i in range(len(value_list)):
        if value_list[i] < 0.0001:
            value_list[i] += 0.01
    color = []
    for name in name_list:
        color.append(group_factor.color_dic[name])
    data = [go.Bar(
        x=name_list,
        y=value_list,
        marker=dict(color=color,
                    line=dict(color=color,
                              width=2)))]

    title = """log10 pvalue of the frequency comparison test showing if %s and %s exons (%s set) 
                are more often regulated by spliceosome factors (%s set)<br>Negative bar: %s exons 
                more often regulated by spliceosome factors<br>Positive bar : %s exons more often regulated 
                by spliceosome factors""" % (pos_reg, neg_reg, group_exon, group_factor_content, neg_reg, pos_reg)

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            title="<-- -log10 pvalue %s area  -:-  log10 pvalue %s area -->" % (neg_reg, pos_reg),
            gridcolor='white',
            gridwidth=1,
            zerolinecolor='rgb(200, 200, 200)',
            zerolinewidth=2,
        ),
        xaxis=dict(title="spliceosome factors"),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        shapes=[dict(type="line", x0=-0.5, y0=math.log10(0.05), x1=len(value_list) - 0.5,
                     y1=math.log10(0.05), line=dict(width=0.5, dash="dash")),
                dict(type="line", x0=-0.5, y0=-math.log10(0.05), x1=len(value_list) - 0.5,
                     y1=-math.log10(0.05), line=dict(width=0.5, dash="dash"))],
        paper_bgcolor='white',
        plot_bgcolor='white',
        showlegend=False
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_VS_%s_exon_%s_spliceosome_%s.html" %
                                      (output, pos_reg, neg_reg, group_exon, group_factor_content),
                                      auto_open=False, validate=False)


def main(filename):
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2",
                                                                 "result/GC_AT_group_regulated_U1_U2/")
    dic_res = file_reader(filename)
    groups_factor = ["all"]
    groups_exon = ["pure"]
    for group_exon in groups_exon:
        for group_factor in groups_factor:
            sub_dic = subsample_result_dic(dic_res, group_exon, group_factor)
            list_values, list_factor_name, pos_reg, neg_reg = value_adapter(sub_dic)
            list_values, list_factor_name = sort_values(list_values, list_factor_name)
            figure_maker(list_values, list_factor_name, pos_reg, neg_reg, group_exon, group_factor, output)


if __name__ == "__main__":
    main(sys.argv[1])
