#!/usr/bin/env python3

# -*- coding: utf-8 -*-



import boxplot_gc_content_maker
import boxplot_flanking_intron_size
import boxplot_gene_size
import os
import sqlite3
import plotly.graph_objs as go
import plotly
import math
import sys
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import numpy as np
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace(
    "boxplot_GC_content_and_flanking_intron_size", ""))
import group_factor


def mann_withney_test_r(list_values1, list_values2):
    """
    Perform a mann withney wilcoxon test on ``list_values1`` and ``list_values2``.

    :param list_values1: (list of float)  list of float
    :param list_values2: (list of float)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
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
    list_values[-1] = list(map(float, list_values[-1]))
    color_dic = group_factor.color_dic
    color_bright = group_factor.color_dic_bright
    data = []
    title = """%s of %s %s containing those exons regulated by different factors""" % (name_fig, regulation, type_fig)
    if "intron" in name_fig:
        list_values = [list(map(math.log10, cur_list)) for cur_list in list_values]
        name_fig = "log " + name_fig
        title += "<br/>(log min flanking intron size for exons and log median introns size for genes)"
    if "gene_size" in name_fig:
        list_values = [list(map(math.log10, cur_list)) for cur_list in list_values]
        name_fig = "log " + name_fig
        title = """%s of %s exons of different exons set""" % (name_fig, regulation)
    for i in range(len(list_name) - 1):
        for j in range(i + 1, len(list_name)):
            pval = mann_withney_test_r(list_values[i], list_values[j])
            title += "<br> mw two-sided test %s vs %s : p = %.2E" % (list_name[i], list_name[j], pval)
    for i in range(len(list_values)):
        cur_color = color_dic["_".join(list_name[i].split("_")[:-1])]
        color_b = color_bright["_".join(list_name[i].split("_")[:-1])]
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True, "fillcolor": color_b, "opacity": 1,
                     "line": {"color": "black"},
                     "box": {"visible": True,  "fillcolor": cur_color}, "meanline": {"visible": False}})

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
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
        showlegend=True,
        shapes=[dict(type="line", x0=-0.5, y0=np.median(list_values[-1]), x1=len(list_values) - 0.5,
                     y1=np.median(list_values[-1]),
                     line=dict(color=color_dic["_".join(list_name[-1].split("_")[:-1])]))]
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_boxplot_%s.html" % (output, name_fig, regulation, type_fig),
                        auto_open=False, validate=False)


def main():
    exon_type = "CCE"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                                 "result/boxplot_gc_content_and_flanking_intron_size/")
    if not os.path.isdir(output):
        os.mkdir(output)
    path = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                               "result/")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/boxplot_GC_content_and_flanking_intron_size",
                                                                "data/sed.db")
    cnx = sqlite3.connect(seddb)
    at_file_pure = "%sAT_rich_exons" % path
    gc_file_pure = "%sGC_rich_exons" % path
    levels = ["exons", "genes"]
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
                    list_gc_content.append(boxplot_gc_content_maker.get_exon_control_gc_content(cnx, exon_type))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.get_exon_control_min_flanking_intron_size(cnx, exon_type)
                    )
            if "genes" in name_file[i]:
                if exon_type not in name_file[i]:
                    list_gc_content.append(boxplot_gc_content_maker.extract_gene_gc_content_from_file(
                        cnx, list_file[i]))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.extract_gene_median_intron_size_from_file(cnx, list_file[i])
                    )
                    list_gene_size.append(boxplot_gene_size.extract_gene_size_from_file(cnx, list_file[i]))
                else:
                    list_gc_content.append(boxplot_gc_content_maker.get_gene_control_gc_content(cnx, exon_type))
                    list_intron_size.append(
                        boxplot_flanking_intron_size.get_gene_control_median_flanking_intron_size(cnx, exon_type)
                    )
                    list_gene_size.append(boxplot_gene_size.get_control_gene_size(cnx, exon_type))
        create_figure(list_gc_content, name_file, output, "down", "GC_content", my_level)
        if my_level == "exons":
            create_figure(list_intron_size, name_file, output, "down", "min_intron_size", my_level)
        if my_level == "genes":
            create_figure(list_intron_size, name_file, output, "down", "median_intron_size", my_level)
            create_figure(list_gene_size, name_file, output, "down", "gene_size", my_level)


if __name__ == "__main__":
    main()
