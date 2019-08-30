#!/usr/bin/env python3.5

"""
Description:
This script contains class and method to compute the GC content of \
exons and the GC content of the other exons contained in the tad.
"""

import numpy as np
import pandas as pd
import os
import plotly
import plotly.graph_objs as go
from scipy import stats
import sys
from functools import reduce


def read_bed(bed_file):
    """
    Get every line of a bed file.

    :param bed_file: (str) a bed file
    :return: (list of list of data) every line in the bed file
    """
    bed_list = []
    with open(bed_file, "r") as bedin:
        for line in bedin:
            if line[0] != "#":
                line = line.replace("\n", "")
                line = line.split("\t")
                line[1] = int(line[1])
                line[2] = int(line[2])
                if len(line) > 6:
                    line[-1] = eval(line[-1])
                bed_list.append(line)
    return bed_list


def intersection_tad(exon, list_tad):
    """
    Get the regions intersecting the current exons.

    :param exon: (list of data) list of data corresponding to an exon \
    in a bed order.
    :param list_tad: (list of list of data) alist of line in a bed file \
    containing tad
    :return: (list of str) list of intersecting region
    """
    exon_start = exon[1]
    exon_stop = exon[2]
    if exon_stop <= exon_start:
        raise ValueError("The size of the exon %s if negative or nul" %
                         exon[3])
    res = []
    for line in list_tad:
        if line[2] > line[1]:
            if  line[0] == exon[0] and \
                    (line[1] <= exon_start <= line[2] and
                     line[1] <= exon_stop <= line[2]):
                res.append(line)
    if len(res) > 1:
        raise NameError("Exon is linked to multiple tad")
    elif len(res) == 0:
        return None
    return res[0]


def intersection_exon(tad, list_exons, exon):
    """
    Get the exons entirely within the current tad.

    :param tad: (list of data) list of data corresponding to a TAD \
    in a bed order.
    :param list_exons: (list of list of data) alist of line in a bed file \
    containing tad
    :param exon: (str ) the exon of interest
    :return: (list of str) list of exon within the current tad
    """
    tad_start = tad[1]
    tad_stop = tad[2]
    if tad_stop <= tad_start:
        raise ValueError("The size of the exon %s if negative or nul" %
                         tad[3])
    res = []
    for line in list_exons:
        if line[2] > line[1]:
            if  line[0] == tad[0] and \
                    (tad_start <= line[1] <= tad_stop and
                     tad_start <= line[2] <= tad_stop and
                    line[3] != exon[3]):
                res.append(line)
    if len(res) == 0:
        return None
    return res


def compute_data(exon_list, tad_list, output, target):
    """
    Return a dataframe showing the ``target`` feature of exons and the \
    average ``target`` feature`` of the other exon in the same tad.

    :param exon_list: (list of list of 7 data) list of exons in bed order
    :param tad_list: (list of list of 6 data) list of tad in bed order
    :param output: (str) folder were the result will be created
    :param target: (str) the target of interest
    :return: (pandas dataframe) table containing the GC content \
    of every exons and the one of other exons in the same tad
    """
    count = 0
    tot = len(exon_list)
    mfile = open("%s/intersection.log" % output, "w")
    ctarget = "%s_other_exons_in_same_tad" % target
    dic = {"exon": [], "tad": [], target: [], ctarget: []}
    for exon in exon_list:
        count += 1
        sys.stdout.write("Processing : %s / %s\t\t\t\r" % (count, tot))
        mfile.write("Exon -> %s\n" % "\t".join(map(str, exon)))
        tad = intersection_tad(exon, tad_list)
        if tad:
            mfile.write("\t|_ TAD -> %s\n" % "\t".join(map(str, tad)))
        else:
            mfile.write("\t|_ TAD -> None\n")
        if tad is not None:
            exons = intersection_exon(tad, exon_list, exon)
            mfile.write("\t\t|_ exons in tad -> %s\n" % str(exons))
            if exons is not None:
                gc_other = np.mean([float(e[-1][target]) for e in exons])
                dic["exon"].append(exon[3])
                dic["tad"].append(tad[3])
                dic[target].append(exon[-1][target])
                dic[ctarget].append(float(gc_other))
    mfile.close()
    df = pd.DataFrame(dic)
    return df[["exon", "tad", target, ctarget]]




def figure_creator_exon(values_xaxis, values_yaxis, exon_name, name_xaxis,
                        name_yaxis, output, name_fig):
    """
    Create a scatter plot showing the potential correlation for every exons \
    regulated each splicing factor in multiple cell lines.

    :param values_xaxis: (list of float) the list of relative value of a \
    particular feature for exons ``regulation`` regulated by ``sf_name``
    :param values_yaxis: (list of float) the list of relative value of the \
    frequency in ``nt`` for exons ``regulation`` regulated by ``sf_name``
    :param exon_name: (list of str) list of exon name
    :param name_xaxis: (string) the name of the features displayed in xaxis
    :param name_yaxis: (string) the name of the features displayed in yaxis
    :param output: (string) path where the results will be created
    :param name_fig: (string) the name of the graphics
    """
    values_xaxis = np.array(values_xaxis, dtype=float)
    values_yaxis = np.array(values_yaxis, dtype=float)
    data = []
    slope, intercept, r_value, p_value, std_err = \
        stats.linregress(values_xaxis, values_yaxis)
    line = slope * np.array(values_xaxis) + intercept
    cor, pval = stats.pearsonr(values_xaxis, values_yaxis)
    data.append(go.Scatter(
        x=values_xaxis,
        y=values_yaxis,
        name="exons",
        mode='markers',
        text=exon_name,
        marker=dict(
            size=7,
            color="navy"
            )
    ))
    data.append(go.Scatter(x=values_xaxis,
                           y=line,
                           mode='lines',
                           line=dict(width=3,
                                     color="red"),
                           name="Fit : p=%.2E, c=%.2f" % (pval, cor)
                           ))

    main_title = 'Correlation between %s and %s ' \
                 '<br> cor : %s - pval : %.2E' \
                 % (name_xaxis, name_yaxis, round(cor, 2), pval)

    x_title = name_xaxis
    y_title = name_yaxis
    figname = '%s/%s.html' % (output, name_fig)
    layout = go.Layout(
        title=main_title,
        hovermode='closest',
        xaxis=dict(
            title=x_title.replace("absolute", ""),
            tickfont=dict(size=25),
            showgrid=False),
        yaxis=dict(
            title=y_title.replace("absolute", ""),
            tickfont=dict(size=25),
            showgrid=False),
        showlegend=True
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=figname, auto_open=False)


def create_figure_for_targets(exon_list, tad_list, output, list_target,
                              name_exon):
    """
    Create on correlation figure for every targets in ``list_target``.

    :param exon_list: (list of list of 7 data) list of exons in bed order
    :param tad_list: (list of list of 6 data) list of tad in bed order
    :param output: (str) folder were the result will be created
    :param list_target: (list of str) list of the targets of interest.
    :param name_exon: (str= the name of the list of exons used
    """
    list_df = []
    for target in list_target:
        print("Computing %s data" % target)
        df = compute_data(exon_list, tad_list, output, target)
        list_df.append(df)
        ctarget = "%s_other_exons_in_same_tad" % target
        print("Creating figure")
        figure_creator_exon(df[ctarget], df[target],
                            df[ctarget],
                            "%s of exons in same tad" % target,
                            "%s exons" % target, output,
                            "correlation_%s_of_%s_VS_exon_same_tad" %
                            (target, name_exon))
    final_df = reduce(lambda left, right:
                      pd.merge(left, right, on=['exon', 'tad']),
                      list_df)
    final_df.to_csv("%s/target_4_%s_and_exons_in_same_TAD.csv" %
                    (output, name_exon), sep="\t", index=False)

def main():
    """
    Create the correlation graphic bewteen the GC content of GC and AT exons \
    and the GC content of other
    """
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.realpath(__file__))))
    output = base + "/result/correlation_GC-AT-exons_TAD/"
    gc_at_exons_file = output + "/GC_content_of_GC-AT_exons.bed"
    tad_file = base + "/data/K562_Lieberman-raw_TADs.hg19.nochr.bed"
    list_target = ["GC_content", "min_flanking_intron_size"]
    gc_at_exon_list = read_bed(gc_at_exons_file)
    tad_list = read_bed(tad_file)
    create_figure_for_targets(gc_at_exon_list, tad_list, output, list_target,
                              "GC-AT_exons")

    cce_exons_file = output + "/GC_content_of_CCE_exons.bed"
    cce_exon_list = read_bed(cce_exons_file)
    create_figure_for_targets(cce_exon_list, tad_list, output, list_target,
                              "GC-AT_exons")

if __name__ == "__main__":
    main()