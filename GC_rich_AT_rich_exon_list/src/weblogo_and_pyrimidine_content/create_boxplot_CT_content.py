#!/usr/bin/python3

# -*- coding utf-8 -*-


import os
import sqlite3
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import sys
import numpy as np
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("/weblogo_and_pyrimidine_content", ""))
import group_factor

# sed_name 2 real_name
sed2real = {"iupac_upstream_intron_adjacent1": "25 nt", "iupac_upstream_intron_adjacent2": "50 nt", "iupac_upstream_intron_proxi": "100 nt"}


def get_control_ct_content(cnx, exon_type, list_targets):
    """
    Get the pyrimidine content of the upstream  of exons ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_type: (string) the type of control exon wanted
    :param list_targets: (list of string) the list columns name where the wanted pyrimidine frequencies is stored
    :return: (dic of list of float) pyrimidine content of the upstream sequence of the control exon set chosen
    """
    cursor = cnx.cursor()
    dic_ct_content = {}
    for target_column in list_targets:
        if exon_type != "ALL":
            query = """SELECT {}
                           FROM sed
                           WHERE exon_type LIKE '%{}%'
                           """.format(target_column, exon_type)
        else:
            query = """SELECT {}
                           FROM sed
                        """.format(target_column)
        cursor.execute(query)
        tuple_list = cursor.fetchall()
        dic_ct_content[sed2real[target_column]] = []
        for iupac in tuple_list:
            if iupac[0] is not None:
                dic_ct_content[sed2real[target_column]].append(iupac[0].split(";")[7])
    # turn tuple into list
    return dic_ct_content



def calculate_ct_content(cnx, gene_id, exon_pos, target_column):
    """
    Get the pyrimidine content of the upstream sequence of exon within the gene identified by ``gene_id`` and \
    having the position ``exon_pos`` in this gene.
    :param cnx: (sqlite3 connect object) connection to sed database
    :param gene_id: (int) the id of the gene containing the interest exons
    :param exon_pos: (int) the position of the exon in the gene ``gene_id``
    :param target_column: (string) the column name where the wanted pyrimidine frequencies is stored
    :return: (float) the pyrimidine content of the upstream sequence of the exon chosen
    """
    cursor = cnx.cursor()
    query = "SELECT {} FROM sed WHERE gene_id = {} and exon_pos = {}".format(target_column, gene_id, exon_pos)
    cursor.execute(query)
    res = cursor.fetchone()
    if res[0] is None:
        return None
    return res[0].split(";")[7]



def extract_ct_content_from_file(cnx, filename, list_targets):
    """
    Get the pyrimidine content of the upstream sequence of 25, 50 and 100 nt of the exon in ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) a file containg exons
    :param list_targets: (string) the list of column names that will be used to extract the y content of upstream \
     sequences of targets exons.
    :return: (dic of list of float) each list corresponds to the pyrimidine content of a particular 25,50,100-nt \
    long upstream sequences of each exon  in ``filename``
    """

    dic_gc = {sed2real[target]:[] for target in list_targets}
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            print(line)
            for i in range(len(list_targets)):
                val = calculate_ct_content(cnx, line[0], line[1], list_targets[i])
                print(val)
                if val is not None:
                    dic_gc[sed2real[list_targets[i]]].append(val)
            line = in_file.readline()
    return dic_gc

def extract_ct_content_from_list(cnx, exon_list, list_targets):
    """

    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_list: (list of 2int) list of exons identified by ``gene_id`` and ``exon_position``
    :param list_targets: (string) the list of column names that will be used to extract the y content of upstream \
     sequences of targets exons.
    :return: (list of list of float) each sublist corresponds to the pyrimidine content of a particular 25,50,100-nt \
    long upstream sequences of each exon  in ``exon_list``
    """
    dic_gc = {sed2real[target]: [] for target in list_targets}
    for exon in exon_list:
        for i in range(len(list_targets)):
            val = calculate_ct_content(cnx, exon[0], exon[1], list_targets[i])
            if val is not None:
                dic_gc[sed2real[list_targets[i]]].append(val)
    return dic_gc


def create_figure(list_values, list_name, output, regulation, name_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    """
    color_dic = group_factor.color_dic
    data = []
    list_values[-1] = list(map(float, list_values[-1]))
    title = """%s of %s exons regulated by different factors""" % (name_fig, regulation)
    for i in range(len(list_values) - 1):
        pval = mann_withney_test_r(list_values[i], list_values[-1])
        title += "<br> %s vs %s, p=%.2E" % (list_name[i], list_name[-1], pval)

    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, # "fillcolor": color_list[i], "opacity": 0.6,
                     "line": {"color": color_dic[list_name[i]]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

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
        shapes=[dict(type="line", x0=-0.5, y0=np.median(list_values[-1]), x1=len(list_values) -0.5, y1=np.median(list_values[-1]),
                     line=dict(color=color_dic[list_name[-1]]))]
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_exons_lvl.html" % (output, name_fig, regulation),
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


def main():

    exon_type = "CCE"
    list_targets = ["iupac_upstream_intron_adjacent1", "iupac_upstream_intron_adjacent2", "iupac_upstream_intron_proxi"]
    path = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "result/")
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "result/boxplot_CT_content/")
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/weblogo_and_pyrimidine_content", "data/sed.db")
    cnx = sqlite3.connect(seddb)

    at_file_pure = "%sAT_rich_exons" % path
    gc_file_pure = "%sGC_rich_exons" % path
    list_file = [gc_file_pure, at_file_pure, None]
    name_file = ["GC_pure", "AT_pure", exon_type]
    list_ct_content = []
    for i in range(len(name_file)):
        if name_file[i] != exon_type:
            list_ct_content.append(extract_ct_content_from_file(cnx, list_file[i], list_targets))
        else:
            list_ct_content.append(get_control_ct_content(cnx, exon_type, list_targets))
    for target in list_targets:
        new_targets_ct_content = []
        for i in range(len(list_ct_content)):
            new_targets_ct_content.append(list_ct_content[i][sed2real[target]])
        create_figure(new_targets_ct_content, name_file, output, "down", "CT content upstream intron %s" %  sed2real[target])



if __name__ == "__main__":
    main()