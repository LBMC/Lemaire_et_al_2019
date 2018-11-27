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


def create_figure(list_values, list_name, output, regulation, name_fig):
    """
    Create a figure showing the values of ``target_name`` for every exon list regulated by ``list_name`` factor

    :param list_values: (list of list of float) each sublist corresponds to the value of  ``target_name`` \
    for every factor in ``list_name``
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param regulation: (string) up or down
    :param name_fig: (string) the name of the graphic
    """
    color_list = ["#0000FF", "#0394d9", "#00aa00", "#03d994", "#FF0000"] * 2
    data = []
    title = """%s of %s exons and gene containing those exons regulated by different factors""" % (name_fig, regulation)
    if "intron" in name_fig:
        list_values = [list(map(math.log10, cur_list)) for cur_list in list_values]
        title += "<br/>(median flanking intron size for exons and median introns size for genes)"
    if "gene_size" in name_fig:
        list_values = [list(map(math.log10, cur_list)) for cur_list in list_values]
        title = """%s of %s exons of different exons set""" % (name_fig, regulation)
    for i in range(len(list_values)):
        data.append({"y": list_values[i], "type": "violin",
                     "name": list_name[i], "visible": True,  # "fillcolor": color_list[i], "opacity": 0.6,
                     "line": {"color": color_list[i]},
                     "box": {"visible": True}, "meanline": {"visible": True}})

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
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_boxplot.html" % (output, name_fig, regulation),
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
    at_file_all = "%sAT_rich_with_intersection_exons" % path
    gc_file_all = "%sGC_rich_with_intersection_exons" % path
    name_file = ["GC_pure_exons", "GC_all_exons", "AT_pure_exons", "AT_all_exons", "%s_exons" % exon_type]
    name_file = name_file + [res.replace("exons", "genes") for res in name_file]
    list_file = [gc_file_pure, gc_file_all, at_file_pure, at_file_all, None] * 2
    list_gc_content = []
    list_intron_size = []
    list_gene_size = []
    for i in range(len(name_file)):
        if "exons" in name_file[i]:
            if exon_type not in name_file[i]:
                list_gc_content.append(boxplot_gc_content_maker.extract_exon_gc_content_from_file(cnx, list_file[i]))
                list_intron_size.append(
                    boxplot_flanking_intron_size.extract_exon_median_flanking_intron_size_from_file(cnx, list_file[i])
                )
            else:
                list_gc_content.append(boxplot_gc_content_maker.get_exon_control_gc_content(cnx, exon_type))
                list_intron_size.append(
                    boxplot_flanking_intron_size.get_exon_control_median_flanking_intron_size(cnx, exon_type)
                )
        if "genes" in name_file[i]:
            if exon_type not in name_file[i]:
                list_gc_content.append(boxplot_gc_content_maker.extract_gene_gc_content_from_file(cnx, list_file[i]))
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
    create_figure(list_gc_content, name_file, output, "down", "GC_content")
    create_figure(list_intron_size, name_file, output, "down", "median_intron_size")
    create_figure(list_gene_size, name_file, output, "down", "gene_size")


if __name__ == "__main__":
    main()
