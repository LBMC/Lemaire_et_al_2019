#!/usr/bin/python3

"""
Description:
    Create a barplot graphics showing if small exons are enriched in GC exons \
    when compared to control exons.
"""

import numpy as np
import plotly.graph_objs as go
import plotly
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from rpy2.robjects.packages import importr
import pandas as pd
import lazyparser as lp


def bed2dic(bed):
    """
    Turn a bed into a dic bed.

    :param bed: (string) a bed file.
    :return: (dictionary of list) dictionary of a bed file containing exons
    """
    dic = {}
    with open(bed, 'r') as bf:
        for line in bf:
            line = line.replace("\n", "").split("\t")
            freq = eval(line[6])
            dic[line[3]] = [line[0], int(line[1]), int(line[2]),
                            line[5], freq]
    return dic



def get_number_of_gc_at_like_exons(dic_exons, ctrl_gc,
                                   ctrl_min_flanking_intron):
    """
    Return the number of GC/AT like exons.

    :param dic_exons: (dictionary of list) dictionary of a bed file \
    containing exons.
    :param ctrl_gc: (float) the median gc content of control exons
    :param ctrl_min_flanking_intron: (float) the min flanking intron size \
     of control exons
    :return: (dictionary) the number of gc like exons at like enxon and total \
    exons in the list.
    """
    dic = {"tot": len(dic_exons.keys()), "gc_like": 0, "at_like": 0}
    for exon in dic_exons.keys():
        gc_content = dic_exons[exon][-1]["GC_content"]
        min_flank_intron = np.nanmin(np.array(
            [dic_exons[exon][-1]["downstream_intron_size"],
             dic_exons[exon][-1]["upstream_intron_size"]], dtype=float))
        if gc_content > ctrl_gc and \
            min_flank_intron < ctrl_min_flanking_intron:
            dic["gc_like"] += 1
        if gc_content < ctrl_gc and \
                min_flank_intron > ctrl_min_flanking_intron:
            dic["at_like"] += 1
    return dic


def get_proportion(dic, typeexons):
    """
    Get the proportion of ``typeexons``.

    :param dic: (dictionary of float)
    :param typefig: (str) the type of exons of interest
    :return: (float)
    """
    return dic[typeexons] / dic["tot"] * 100


def create_barplot(list_proportion, list_name, output, name_fig, typefig):
    default_colors = ["green", "blue", "purple"]
    data = []
    lcol = []
    for i in range(len(list_proportion)):
        if "CCE" in list_name[i]:
            mcol = "red"
        else:
            mcol = default_colors[i]
        lcol.append(mcol)
    data.append(go.Bar(x=list_name, y=list_proportion,
                       name=list_name, marker=dict(color=lcol)))

    title = """Proportion of %s exons""" % typefig

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
    plotly.offline.plot(fig, filename="%s/%s_for_%s_exons.html" %
                                      (output, name_fig, typefig),
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

    rm1 = tot1 - obs1
    rm2 = tot2 - obs2
    vect = v.FloatVector([obs1, rm1, obs2, rm2])
    pval = float(chisq(vect)[0])
    print(obs1, rm1, obs2, rm2)
    print(pval)
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


def write_proportion_pvalues(list_dic_prop, list_name, output, name_fig):
    """
    Write a text file containing the pvalue of a frequency test for test sets vs control. \

    :param list_dic_prop: (list of dic of float) the number of gc like exons at like enxon and total \
    exons in the list.
    :param list_name: (list of string) list of name of different factor studied
    :param output: (string) path where the output_file will be created
    :param name_fig: (string) the name of the figure
    """
    list_col = ["Categorie", "factor1", "ctrl", "nb_like1", "tot1",
                "nb_like_ctrl",  "tot_ctrl", "pval"]
    data = {my_name: [] for my_name in list_col}
    filename = "%s/%s.txt" % (output, name_fig)
    for i in range(len(list_dic_prop) - 1):
        for typefig in ["gc_like", "at_like"]:
            tot1 = list_dic_prop[i]["tot"]
            tot2 = list_dic_prop[-1]["tot"]
            nb1 = list_dic_prop[i][typefig]
            nb2 = list_dic_prop[-1][typefig]
            pval = frequency_test(nb1, tot1, nb2, tot2)
            list_info = [typefig, list_name[i], list_name[-1], nb1,
                         tot1, nb2, tot2, pval]
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


@lp.parse(list_bed="file", output="dir")
def main(list_bed, bed_names, output):
    """
    Create two barplots figure to see if small exons are enriched
    in GC-like exons.

    :param list_bed: (List(vtype=str)) list of bed file
    :param bed_names: (List(vtype=str)) list of bed names
    :param output: (str) path were the result will be created
    """
    ctrl_gc = 49.3
    ctrl_mfi = 691
    list_dic_prop = [get_number_of_gc_at_like_exons(bed2dic(mfile),
                                                    ctrl_gc, ctrl_mfi)
                     for mfile in list_bed]
    list_prop_gc = [get_proportion(mdic, "gc_like") for mdic in list_dic_prop]
    list_prop_at = [get_proportion(mdic, "at_like") for mdic in list_dic_prop]
    name_fig = "barplot_small_big_exons"
    create_barplot(list_prop_at, bed_names, output, name_fig, "at_like")
    create_barplot(list_prop_gc, bed_names, output, name_fig, "gc_like")
    write_proportion_pvalues(list_dic_prop, bed_names, output, name_fig)


if __name__ == "__main__":
    main()


