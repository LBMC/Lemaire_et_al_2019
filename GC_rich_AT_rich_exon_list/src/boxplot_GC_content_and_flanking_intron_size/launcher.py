#!/usr/bin/env python3

# -*- coding: utf-8 -*-



import boxplot_gc_content_maker
import boxplot_flanking_intron_size
import boxplot_gene_size
import os
import sqlite3
import copy
import plotly.graph_objs as go
import plotly
import math
import sys
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stat_maker
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace(
    "boxplot_GC_content_and_flanking_intron_size", ""))
import group_factor
import union_dataset_function



def mann_withney_test_r(list_values1, list_values2, alt="less"):
    """
    Perform a mann withney wilcoxon test on ``list_values1`` and ``list_values2``.

    :param list_values1: (list of float)  list of float
    :param list_values2: (list of float)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative="%s", correct=F)
        return(test$p.value)
    }

                   """ % alt)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def create_figure(list_values, list_name, output, regulation, name_fig, type_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param name_fig: (string) the name of the graphic
    :param type_fig: (string) the type of graphic build
    """
    list_values2 = copy.deepcopy(list_values)
    for i in range(len(list_values2)):
        v = np.array(list_values2[i], dtype=float)
        list_values2[i] = list(v[~np.isnan(v)])
    color_dic = group_factor.color_dic
    color_bright = group_factor.color_dic_bright
    data = []
    title = """%s of %s %s containing those exons regulated by different factors""" % (name_fig, regulation, type_fig)
    if "intron" in name_fig:
        list_values2 = [list(map(math.log10, cur_list)) for cur_list in list_values2]
        name_fig = "log " + name_fig
        title += "<br/>(log min flanking intron size for exons and log median introns size for genes)"
    if "gene_size" in name_fig:
        list_values2 = [list(map(math.log10, cur_list)) for cur_list in list_values2]
        name_fig = "log " + name_fig
        title = """%s of %s exons of different exons set""" % (name_fig, regulation)
    if type_fig == "exons" or name_fig == "log median_intron_size":
        for i in range(len(list_name) - 1):
            for j in range(i + 1, len(list_name)):
                pval1 = mann_withney_test_r(list_values2[i], list_values2[j], alt="two.sided")
                title += "<br> mw 2-sided test %s vs %s : p = %.2E" % (list_name[i], list_name[j], pval1)
    for i in range(len(list_values2)):
        cur_color = color_dic["_".join(list_name[i].split("_")[:-1])]
        color_b = color_bright["_".join(list_name[i].split("_")[:-1])]
        if type_fig == "exons":
            data.append({"y": list_values2[i], "type": "violin",
                         "name": list_name[i], "visible": True, "fillcolor": color_b, "opacity": 1,
                         "line": {"color": "black"},
                         "box": {"visible": True,  "fillcolor": cur_color}, "meanline": {"visible": False}})
        else:
            data.append({"y": list_values2[i], "type": "box",
                         "name": list_name[i], "visible": True, "fillcolor":cur_color , "opacity": 1,
                         "line": {"color": "black"}})

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            title=name_fig,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
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
        showlegend=True
        # shapes=[dict(type="line", x0=-0.5, y0=np.median(list_values[-1]), x1=len(list_values) - 0.5,
        #              y1=np.median(list_values[-1]),
        #              line=dict(color=color_dic["_".join(list_name[-1].split("_")[:-1])]))]
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_boxplot_%s.html" % (output, name_fig, regulation, type_fig),
                        auto_open=False, validate=False)


def dataframe_creator(list_values, list_name, output, regulation, name_df, type_fig):
    """
    Create a pandas dataframe.

    :param list_values: (list of list of float) list of values
    :param list_name: (list of string) the names of the exon list used to get each sublist on ``list_values`` objects.
    :param name_df: (string) the type of values displayed in ``list_values``
    :param output: (string) folder where the result will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of graphic build
    :return: (pandas Dataframe) a dataFrame with the values in ``list_values`` and the names in ``list_names``
    """
    new_values = np.hstack(list_values)
    new_names = np.hstack([[list_name[i]] * len(list_values[i]) for i in range(len(list_values))])
    new_values = new_values.astype(np.float)
    new_names = list(new_names[~np.isnan(new_values)])
    new_values = list(new_values[~np.isnan(new_values)])
    df = pd.DataFrame({"values": new_values, "project": new_names})
    filename = "%s%s_%s_dataframe_%s_stat.txt" % (output, name_df, regulation, type_fig)
    # if name_df == "GC_content" and type_fig == "genes":
    #     new_df = stat_maker.nb_glm_stats(df, filename)
    #     new_df.to_csv("%s%s_%s_dataframe_%s_stat.txt" % (output, name_df, regulation, type_fig), index=False, sep="\t")
    if name_df == "gene_size":
        df["values"] = np.log10(df["values"].values)
        new_df = stat_maker.anova_nt_stats(df, filename)
        new_df.to_csv("%s%s_%s_dataframe_%s_stat.txt" % (output, name_df, regulation, type_fig), index=False, sep="\t")
    df.to_csv("%s%s_%s_dataframe_%s_table.txt" % (output, name_df, regulation, type_fig), index=False, sep="\t")



def dataframe_creator2(list_values, list_values2, list_name, output, regulation, name_df, name_df2, type_fig):
    """
    Create a pandas dataframe.

    :param list_values: (list of list of float) list of values
    :param list_name: (list of string) the names of the exon list used to get each sublist on ``list_values`` objects.
    :param name_df: (string) the type of values displayed in ``list_values``
    :param output: (string) folder where the result will be created
    :param regulation: (string) up or down
    :param type_fig: (string) the type of graphic build
    :return: (pandas Dataframe) a dataFrame with the values in ``list_values`` and the names in ``list_names``
    """
    new_values = np.hstack(list_values)
    new_values2 = np.hstack(list_values2)
    new_names = np.hstack([[list_name[i]] * len(list_values[i]) for i in range(len(list_values))])
    new_values = new_values.astype(np.float)
    new_values2 = new_values2.astype(np.float)
    df = pd.DataFrame({name_df: new_values, "project": new_names, name_df2: new_values2})
    filename = "%s%s-%s_%s_dataframe_%s_table.txt" % (output, name_df, name_df2, regulation, type_fig)
    df.to_csv(filename, index=False, sep="\t")
    new_df = stat_maker.anova_gene_stats(df, filename, name_df, name_df2)
    new_df.to_csv(filename.replace("table.txt", "stat.txt"), index=False, sep="\t")


def extract_gene(filename):
    """
    From a file with two columns, the id of a gene and the exon number, return
    :param filename: (string ) a file of exons
    :return: (a list of string) the list of gene id
    """
    genes = []
    with open(filename, "r") as my_file:
        for line in my_file:
            line = line.split("\t")[0]
            genes.append(line)
    genes = list(np.unique(genes))
    return genes


def venn_diagram_creator(list1, name1, list2, name2, output):
    """
    Create a venn diagram of the `list1` and `list2` lists of values.

    :param list1: (list of 2 int) list of gene identified by their gene_id
    :param name1: (string) the name of the list of gene named `list1`
    :param list2: (list of 2 int) list of gene identified by their gene_id
    :param name2: (string) the name of the list of gene named `list2`
    :param output: (string) path where the result venn diagram will be created
    :return:
    """
    # plt.figure(figsize=(48. / 2.54, 27 / 2.54))
    venn2([set(list1), set(list2)], set_labels = (name1, name2))
    plt.savefig("%sVenn_%s_vs_%s.pdf" % (output, name1, name2))
    plt.clf()
    plt.cla()
    plt.close()


def get_common_genes(filename1, filename2, output):
    """
    From 2 files containing list of exons : return the common genes between those 2 lists
    :param filename1: (string) a file containing a list of exons
    :param filename2: (string) a file containing a list of exons
    :return: (list of string) the list genes common between ``filename1``, ``filename2``
    """
    genes1 = extract_gene(filename1)
    genes2 = extract_gene(filename2)
    venn_diagram_creator(genes1, "AT_genes", genes2, "GC_genes", output)
    return [g for g in genes1 if g in genes2]



def main():
    exon_type = "CCE"
    regulation = "down"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                                 "result/boxplot_gc_content_and_flanking_intron_size/")
    if not os.path.isdir(output):
        os.mkdir(output)
    path = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                               "result/")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                                "data/sed.db")
    cnx = sqlite3.connect(seddb)
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx, regulation)
    gene2remove = [exon[0] for exon in exon2remove]
    at_file_pure = "%sAT_rich_exons" % path
    gc_file_pure = "%sGC_rich_exons" % path
    gene2remove_at_gc = get_common_genes(at_file_pure, gc_file_pure, output)
    #levels = ["exons", "genes"]
    levels = ["genes"]
    for my_level in levels:
        name_file = ["GC_pure_%s" % my_level, "AT_pure_%s" % my_level, "%s_%s" % (exon_type, my_level)]
        list_file = [gc_file_pure, at_file_pure, None]
        list_gc_content = []
        list_intron_size = []
        list_gene_size = []
        for i in range(len(name_file)):
            if "exons" in name_file[i]:
                if exon_type not in name_file[i]:
                    list_gc_content.append(boxplot_gc_content_maker.extract_exon_gc_content_from_file(
                        cnx, list_file[i]))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.extract_exon_min_flanking_intron_size_from_file(cnx, list_file[i])
                    )
                else:
                    list_gc_content.append(boxplot_gc_content_maker.get_exon_control_gc_content(cnx, exon_type, exon2remove))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.get_exon_control_min_flanking_intron_size(cnx, exon_type, exon2remove)
                    )
            if "genes" in name_file[i]:
                print(name_file[i])
                if exon_type not in name_file[i]:
                    list_gc_content.append(boxplot_gc_content_maker.extract_gene_gc_content_from_file(
                        cnx, list_file[i], gene2remove_at_gc))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.extract_gene_median_intron_size_from_file(cnx, list_file[i], gene2remove_at_gc)
                    )
                    list_gene_size.append(boxplot_gene_size.extract_gene_size_from_file(cnx, list_file[i], gene2remove_at_gc))
                else:
                    list_gc_content.append(boxplot_gc_content_maker.get_gene_control_gc_content(cnx, exon_type, gene2remove))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.get_gene_control_median_flanking_intron_size(cnx, exon_type, gene2remove)
                    )
                    list_gene_size.append(boxplot_gene_size.get_control_gene_size(cnx, exon_type, gene2remove))
        if my_level == "exons":
            create_figure(list_gc_content, name_file, output, regulation, "GC_content", my_level)
            dataframe_creator(list_gc_content, name_file, output, regulation, "GC_content", my_level)
            create_figure(list_intron_size, name_file, output, regulation, "min_intron_size", my_level)
            dataframe_creator(list_intron_size, name_file, output, regulation, "min_intron_size", my_level)
        if my_level == "genes":
            # create_figure(list_intron_size, name_file, output, regulation, "median_intron_size", my_level)
            # dataframe_creator(list_intron_size, name_file, output, regulation, "median_intron_size", my_level)
            # create_figure(list_gene_size, name_file, output, regulation, "gene_size", my_level)
            # dataframe_creator(list_gene_size, name_file, output, regulation, "gene_size", my_level)
            # create_figure(list_gc_content, name_file, output, regulation, "GC_content", my_level)
            dataframe_creator2(list_gc_content, list_gene_size, name_file, output, regulation, "GC_content", "gene_size", my_level)
            dataframe_creator2(list_intron_size, list_gc_content, name_file, output, regulation,
                               "median_intron_size", "GC_content", my_level)


if __name__ == "__main__":
    main()
