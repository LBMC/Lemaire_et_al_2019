#!/usr/bin/env python3.5

# -*- coding: UTF-8 -*-

"""
Description:
    The goal of this script is to determine whether, in average, the TADs \
    contains exons with an homogeneous number of branch point in their \
    upstream intron or an heterogeneous number.
"""

import pandas as pd
import tad_exon as te
import lazyparser as lp
import math
import numpy as np
import seaborn as sns
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()


def intersection_exon(tad, list_exons):
    """
    Get the exons entirely within the current tad.

    :param tad: (list of data) list of data corresponding to a TAD \
    in a bed order.
    :param list_exons: (list of list of data) alist of line in a bed file \
    containing tad
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
                     tad_start <= line[2] <= tad_stop):
                res.append(line)
    if len(res) == 0:
        return None
    return res


def filter_dataframe(dataframe, threshold=10):
    """
    Remove every tad that doesn't contain at least ``threshold`` exons.

    :param dataframe: (pandas dataframe) a dataframe of tad
    :param threshold: (int) a number
    :return: (pandas dataframe) the dataframe of tad having at least \
    ``threshold`` exons
    """
    tads, counts = np.unique(dataframe["Tad"], return_counts=True)
    tad_2_keep = []
    for tad, count in zip(tads, counts):
        if count >= threshold:
            tad_2_keep.append(tad)
    dataframe = dataframe[(dataframe["Tad"].isin(tad_2_keep))]
    return dataframe


def get_other_regulated(exon, list_exons, target):
    """
    Get the number of exons regulated by ``target`` in ``list_exons`` \
    after removing ``exon`` from that list.

    :param exon: (list of values) an exon
    :param list_exons: (list of list of value) a list of exons
    :return: (int) the number of exons in ``list_exons`` regulated \
    by target
    """
    start = len(list_exons)
    list_exons = [mexon for mexon in list_exons if mexon != exon]
    stop = len(list_exons)
    if start - stop != 1:
        raise ValueError("Exon %s must be in %s" % (exon, list_exons))
    return sum([mexon[-1][target] for mexon in list_exons])


def compute_dataframe(exon_list, tad_list, output, name_table):
    """
    Return a dataframe showing the exons within a tad and their content \
    in branch points in their upstream intron.

    :param exon_list: (list of list of 7 data) list of exons in bed order
    :param tad_list: (list of list of 6 data) list of tad in bed order
    :param output: (str) folder were the result will be created
    :param name_table: (str) the name of the table to create
    :return: (pandas dataframe) table containing the tad \
    the exons it contains and their number of branch point in  \
    their upstream intron
    """
    count = 0
    dic = {"Tad": [], "Exon": [], "good_bp": [], "bp>2": [], "size_tad": [],
           "U1_regulated": [], "U2_regulated": [], "UNA_count": [],
           "T_stretch": [], "MFE_5SS": [], "MFE_3SS": [], "other_U1": [],
           "other_U2": []}
    for tad in tad_list:
        exons = intersection_exon(tad, exon_list)
        if exons is None:
            count += 1
        else:
            for exon in exons:
                ou1 = get_other_regulated(exon, exons, "U1-regulated")
                ou2 = get_other_regulated(exon, exons, "U2-regulated")
                cbp = exon[-1]["good_bp"]
                if cbp is not None:
                    dic["Tad"].append(tad[3])
                    dic["Exon"].append(exon[3])
                    dic["good_bp"].append(cbp)
                    dic["bp>2"].append(1 if cbp > 2 else 0)
                    dic["size_tad"].append(len(exons))
                    dic["U1_regulated"].append(exon[-1]["U1-regulated"])
                    dic["U2_regulated"].append(exon[-1]["U2-regulated"])
                    dic["UNA_count"].append(exon[-1]["UNA_count"])
                    dic["T_stretch"].append(exon[-1]["T_stretch"])
                    dic["MFE_5SS"].append(math.sqrt(
                        (exon[-1]["MFE_5SS"] * -1) + 3/8))
                    dic["MFE_3SS"].append(math.sqrt(
                        (exon[-1]["MFE_3SS"] * -1) + 3/8))
                    dic["other_U1"].append(ou1)
                    dic["other_U2"].append(ou2)
    df = pd.DataFrame(dic)
    df = df[["Tad", "Exon", "good_bp", "bp>2", "size_tad", "U1_regulated",
             "U2_regulated", "UNA_count", "T_stretch", "MFE_5SS", "MFE_3SS",
             "other_U1", "other_U2"]]
    df.to_csv("%s/%s.txt" % (output, name_table), sep="\t", index=False)
    print("%s / %s Tads don't contain any (usable) exons" %
          (count, len(tad_list)))
    # df_filtered = filter_dataframe(df, threshold)
    # filename = "%s/%s_with_more_than_%s_exons.txt" % (output, name_table,
    #                                                   threshold)
    # df_filtered.to_csv(filename, sep="\t", index=False)
    return df, "%s/%s.txt" % (output, name_table)


def glmm_analysis(dataframe, filename, output, target, family="binomial"):
    """
    Perform a glmm analysis of the data of interest.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string) list of float
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require(lme4)
    function(data, name, target){
        mod <- glmer(%s ~ size_tad + (1|Tad), family=%s, data=data, nAGQ=0)
        nulmod <- glm(%s ~ size_tad, family=%s, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mod, n = 250)
        png(paste(name, "/mod_dignostics_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        simulationOutput <- simulateResiduals(fittedModel = nulmod, n = 250)
        png(paste(name, "/nulmod_dignostics_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        return(anova(mod, nulmod, test="Chisq"))
    }

    """ % (target, family, target, family, target, target))

    res = stat_s(r_df, output, target)
    print(res)
    pandas2ri.ri2py(res).to_csv("%s_glmm_stats.txt" % filename,
                                sep="\t", index=False)


def glmmnb_analysis(dataframe, filename, output, target):
    """
    Perform a glmm analysis of the data of interest.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string) list of float
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require(lme4)
    require("MASS")
    function(data, name, target){
        mod <- glmer.nb(%s ~ size_tad + (1|Tad), data=data, nAGQ = 0)
        nulmod <- glm.nb(%s ~ size_tad, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mod, n = 250)
        png(paste(name, "/mod_dignostics_%s.png", sep=""), 
            height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        simulationOutput <- simulateResiduals(fittedModel = nulmod, n = 250)
        png(paste(name, "/nulmod_dignostics_%s.png", sep=""), 
            height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        return(anova(mod, nulmod, test="Chisq"))
    }

    """ % (target, target, target, target))

    res = stat_s(r_df, output, target)
    print(res)
    pandas2ri.ri2py(res).to_csv("%s_glmm_stats.txt" % filename,
                                sep="\t", index=False)


def glmm_spliceosome(dataframe, filename, output, target):
    """
    Perform a glmm analysis of the data of interest.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string) list of float
    """
    if target == "U1_regulated":
        variable = "other_U1"
    elif target == "U2_regulated":
        variable = "other_U2"
    else:
        raise ValueError("target must be in [U1_regulated, U2_regulated]"
                         "not %s" % target)
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require(lme4)
    function(data, name, target){
        mod <- glm(%s ~ size_tad + %s, family=binomial, data=data)
        nulmod <- glm(%s ~ size_tad, family=binomial, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mod, n = 250)
        png(paste(name, "/mod_dignostics_spl_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        simulationOutput <- simulateResiduals(fittedModel = nulmod, n = 250)
        png(paste(name, "/nulmod_dignostics_spl_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        return(anova(nulmod, mod, test="Chisq"))
    }

    """ % (target, variable, target, target, target))

    res = stat_s(r_df, output, target)
    print(res)
    pandas2ri.ri2py(res).to_csv("%s_glmm_spliceosome_stats.txt" % filename,
                                sep="\t", index=False)


def lmm_analysis(dataframe, filename, output, target):
    """
    Perform a glmm analysis of the data of interest.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string) list of float
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require(lme4)
    require("MASS")
    function(data, name, target){
        mod <- lmer(%s ~ size_tad + (1|Tad), data=data)
        nulmod <- lm(%s ~ size_tad, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mod, n = 250)
        png(paste(name, "/mod_dignostics_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        simulationOutput <- simulateResiduals(fittedModel = nulmod, n = 250)
        png(paste(name, "/nulmod_dignostics_%s.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        return(anova(mod, nulmod, test="Chisq"))
    }

    """ % (target, target, target, target))

    res = stat_s(r_df, output, target)
    print(res)
    pandas2ri.ri2py(res).to_csv("%s_glmm_stats.txt" % filename,
                                sep="\t", index=False)


def create_figure(df, output, target, col):
    """
    Create a line plot.

    :param df: (pandas Dataframe) a table of value
    :param output: (str) path where the figure will be created
    :param target: (str) the column of interest
    :param col: (list of str) list of color
    """
    df = df[df["Percent"] != "[0:20["]
    df = df.sort_values("Percent")
    sns.set()
    # sns.set_context("poster")
    g = sns.catplot(x="Percent", y="Tad", hue="Name", data=df,
                    height=9, aspect=1.7, ci="sd", palette=col, kind="bar")
    g.savefig("%s/%s.pdf" % (output, target))


def my_apply(x):
    # step = 20
    # intervals = [[i, i+ step] for i in range(0, 100, step)]
    intervals =[[0, 20], [20, 100]]
    intervals[-1][-1] += 1
    for interval in intervals:
        if interval[0] <= x < interval[1]:
            return "[" + str(interval[0]) + ":" + str(interval[1]) + "["
    raise ValueError("x must be between 0 and 100")


def figure_regulated(df, target, output, col, iteration=1000):
    """
    Create the line plot of U1-regulated exons or U2-regulated exons \
    compared to ``iteration`` shuffles.

    :param df: (pandas Dataframe) a pandas dataframe
    :param target: (str) the column of interest
    :param output: (str) path where the figure will be created
    :param col: (list of str) list of color
    :param iteration: (int) the number of shuffles to do
    """
    full_df = pd.DataFrame()
    col_shuf = target.replace("regulated", "shuffled")
    target_percent = target + "_percent"
    col_shuf_percent = col_shuf + "_percent"
    for i in range(iteration):
        df[col_shuf] = df[target].sample(frac=1).reset_index(drop=True)
        dfs = df[["Tad", "size_tad"]].drop_duplicates()
        dfshuf = df.groupby(["Tad"]).sum().reset_index()[
            ["Tad", col_shuf, target]]
        dfr = pd.merge(dfs, dfshuf, on="Tad")
        dfr[target_percent] = round(dfr[target] / dfr["size_tad"] * 100, 0)
        dfr[col_shuf_percent] = round(dfr[col_shuf] / dfr["size_tad"] * 100, 0)

        if i == 0:
            tmp = dfr[[target_percent, "Tad"]].groupby(target_percent).count() / dfr.shape[0] * 100
            tmp["Name"] = [target] * tmp.shape[0]
            tmp["Percent"] = tmp.index
            tmp = tmp.reset_index(drop=True)
            tmp["Percent"] = tmp["Percent"].apply(my_apply)
            tmp = tmp.groupby(["Percent", "Name"]).sum().reset_index()
            print(tmp)
            full_df = pd.concat([full_df, tmp], axis=1)
        tmp = dfr[[col_shuf_percent, "Tad"]].groupby(col_shuf_percent).count() / \
              dfr.shape[
                  0] * 100
        tmp["Name"] = [col_shuf] * tmp.shape[0]
        tmp["Percent"] = tmp.index
        tmp = tmp.reset_index(drop=True)
        tmp["Percent"] = tmp["Percent"].apply(my_apply)
        tmp = tmp.groupby(["Percent", "Name"]).sum().reset_index()
        full_df = pd.concat([full_df, tmp], axis=0)
    create_figure(full_df, output, target, col)


@lp.parse(tad_bed="file", exon_bed="file", output="dir")
def main(tad_bed, exon_bed, output, name_table, threshold=10):
    """
    Determine whether, in average, the TADs contains exons with an homogeneous
    number of branch point in their upstream intron or an heterogeneous number.

    :param tad_bed: (str) A bed file containing Tad
    :param exon_bed: (str) A bed file containing exons and an \
    additional column, corresponding to a python dictionary, with the key \
    'goob_bp'
    """
    exon_list = te.read_bed(exon_bed)
    tad_list = te.read_bed(tad_bed)
    df_filtered, filename = \
        compute_dataframe(exon_list, tad_list, output, name_table)
    glmm_analysis(df_filtered, filename.replace(".txt", "_target_bp>2"), output, "bp.2")
    glmm_analysis(df_filtered, filename.replace(".txt", "_target_U1-reg"),
                  output, "U1_regulated")
    glmm_analysis(df_filtered, filename.replace(".txt", "_target_U2-reg"),
                  output, "U2_regulated")
    # # glmm_analysis(df_filtered, filename.replace(".txt", "_target_UNA"),
    # #               output, "UNA_count", "poisson")
    glmmnb_analysis(df_filtered, filename.replace(".txt", "_target_UNA"),
                  output, "UNA_count")
    glmmnb_analysis(df_filtered, filename.replace(".txt", "_target_T-stretch"),
                  output, "T_stretch")
    # # glmm_analysis(df_filtered, filename.replace(".txt", "_target_T-stretch"),
    # #               output, "T_stretch", "poisson")
    lmm_analysis(df_filtered, filename.replace(".txt", "_MFE_5SS"),
                  output, "MFE_5SS")
    lmm_analysis(df_filtered, filename.replace(".txt", "_MFE_3SS"),
                  output, "MFE_3SS")
    glmm_spliceosome(df_filtered, filename.replace(".txt", "_target_U1-reg"),
                     output, "U1_regulated")
    glmm_spliceosome(df_filtered, filename.replace(".txt", "_target_U2-reg"),
                     output, "U2_regulated")
    figure_regulated(df_filtered, "U2_regulated", output, ["#66FF66", "red"],
                     iteration=100)
    figure_regulated(df_filtered, "U1_regulated", output, ["#9999FF", "red"],
                     iteration=100)

if __name__ == "__main__":
    main()