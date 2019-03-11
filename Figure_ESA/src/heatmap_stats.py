#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import numpy as np
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
import rpy2.robjects.vectors as v
import pandas as pd
pandas2ri.activate()


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


def perform_mann_withney_test(dataframe, sf_name, exon_type):
    """
    From a dataframe of value perform a Mann Withney Wilcoxon test.

    :param dataframe: (pandas DataFrame)
    :param sf_name: (string)
    :param exon_type: (string)
    :return: (pandas dataFrame)
    """
    rstats = importr('stats')
    list_ctrl = np.array(dataframe[dataframe["project"] == exon_type].loc[:, "values"].values, dtype=float)
    list_ctrl = list(list_ctrl[~np.isnan(list_ctrl)])
    pval_list = []
    for my_sf in sf_name:
        test_list = np.array(dataframe[dataframe["project"] == my_sf].loc[:, "values"].values, dtype=float)
        test_list = list(test_list[~np.isnan(test_list)])
        pval_list.append(mann_withney_test_r(test_list, list_ctrl))
    pcor = rstats.p_adjust(v.FloatVector(pval_list), method="BH")
    df = pd.DataFrame({"SF": sf_name, "pval_MW": pval_list, "p_adj": pcor})
    return df[["SF", "pval_MW", "p_adj"]]


def handle_dataframe_statistics(dataframe, filename):
    """
    Make analysis on the dataFrame ``dataframe``

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string) the name of the output file
    :return: (pandas DataFrame) the statistical analyzes
    """
    features = np.unique(dataframe["features"].values)
    val = [n for n in features if "S_nt" in n]
    val2 = [n for n in features if "nt" in n]
    adj = [n for n in features if "adjacent1" in n]
    intron = [n for n in features if "intron" in n]
    if len(features) == 2 and len(val) == 1 and len(intron) == 0:
        dataframe = dataframe[dataframe["features"].str.contains("S_nt")]
        dataframe = S_nt_stats(dataframe, filename)
        dataframe["project"] = [val.replace("project", "") for val in dataframe.index.values]
        dataframe = dataframe[dataframe["project"] != "(Intercept)"] # removing the intercept line
        dataframe = dataframe[["project", "Pr(>|t|)"]]
    if len(features) == 2 and len(val) == 1 and len(intron) > 1:
        dataframe = dataframe[dataframe["features"].str.contains("S_nt")]
        exon_type =  [et for et in np.unique(dataframe["project"].values) if et in ["ACE", "CCE", "ALL"]][0]
        sf_name =  [et for et in np.unique(dataframe["project"].values) if et != exon_type]
        dataframe = perform_mann_withney_test(dataframe, sf_name, exon_type)
    if len(features) == 4 and len(val2) == 4 and len(adj) == 4:
        cor_feature = {f: f.split("_")[0] for f in features}
        dataframe["features"] = dataframe["features"].map(cor_feature)
        dataframe.to_csv(filename.split(".")[0] + "_table.txt", sep="\t", index=False)
        dataframe = nt_count_stats(dataframe, filename)
    return dataframe


def nt_count_stats(dataframe, filename):
    """
    Perform a mann withney wilcoxon test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    function(data, name){
        data$values <- round((data$values / 100) * 25)
        factors <- levels(as.factor(as.vector(data$features)))
        png(paste(name, "_hist.png", sep=""), height=2160, width=1920)
        par(mfrow = c(2, 2))
        for (f in factors){
            hist(data$values[data$features == f],  breaks=sqrt(length(data$values[data$features == f])), main=paste("hist for", f))
        }
        dev.off()
        final <- NULL
        for(f in factors){
            mglm <- glm(values ~ project, data = data[data$features == f, ], family="poisson")
            simulationOutput1 <- simulateResiduals(fittedModel = mglm, n = 250)
            png(paste(name, "_hist_", f, ".png", sep=""), height=2160, width=1920)
            plot(simulationOutput1)
            dev.off()
            tmp <- summary(mglm)$coefficients
            row.names(tmp) <- paste(row.names(tmp), ".", f, sep="")
            final <- rbind(final, tmp)
        }
        return(as.data.frame(final))
    }

                   """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    df_stats["index"] = df_stats.index
    return df_stats


def S_nt_stats(dataframe, filename):
    """
    Perform a mann withney wilcoxon test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""

    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        png(paste(name, "_qqplots.png", sep=""), height=2160, width=1920)
        par(mfrow = c(5, 7))
        pvalues <- NULL
        i = 0
        for (f in factors){
            i = i + 1
            test <- ks.test(data$values[data$project == f], "pnorm", mean=mean(data$values[data$project == f]), 
                            sd=sd(data$values[data$project == f]))
            pvalues[i] <- test$p.value
            qqnorm(data$values[data$project == f], main=paste("normal QQ-PLOT for", f, "ks-test :", test$p.value))
        }
        dev.off()

        mlm <- lm(values ~ project, data=data)
        png(paste(name, "_dignostics.png", sep=""), height=2160, width=870)
        par(mfrow=c(4, 1))
        plot(mlm)
        dev.off()

        return(as.data.frame(summary(mlm)$coefficients))
    }

                   """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    return df_stats