#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import numpy as np
from rpy2.robjects import r, pandas2ri

pandas2ri.activate()


def anova_nt_stats(dataframe, filename):
    """
    Perform an anova test on ``dataframe``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""

    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        png(paste(name, "_distrib_qqplots.png", sep=""), height=2160, width=1920)
        par(mfrow = c(3, 2))
        pvalues <- NULL
        i = 0
        for (f in factors){
            i = i + 1
            test <- ks.test(data$values[data$project == f], "pnorm", mean=mean(data$values[data$project == f]),
                            sd=sd(data$values[data$project == f]))
            pvalues[i] <- test$p.value
            qqnorm(data$values[data$project == f], main=paste("normal QQ-PLOT for", f, "ks-test :", test$p.value))
            hist(data$values[data$project == f],  breaks = sqrt(length(data$values[data$project == f])),
                 main=paste("hist for", f))
        }
        dev.off()

        mlm <- aov(values ~ project, data=data)
        png(paste(name, "_dignostics.png", sep=""), height=2160, width=870)
        par(mfrow=c(2, 2))
        plot(mlm)
        dev.off()
        return(as.data.frame(TukeyHSD(mlm)$project))
    }

                   """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    df_stats["project"] = df_stats.index.values
    df_stats = df_stats[["project", "diff", "lwr", "upr", "p adj"]]
    return df_stats


def nb_glm_stats(dataframe, filename):
    """
    Perform a glm nb test on ``dataframe``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    dataframe = dataframe[dataframe["project"] != "CCE_exons"]
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("MASS")
    require("DHARMa")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        png(paste(name, "_distrib.png", sep=""), height=2160, width=1920)
        par(mfrow = c(2, 2))
        for (f in factors){
            hist(log10(data$values[data$project == f]),  breaks = sqrt(length(data$values[data$project == f])),
            main=paste("hist for", f))
        }
        dev.off()
        mlm <- glm.nb(values ~ project, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mlm, n = 250)
        png(paste(name, "_dignostics.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()

        return(as.data.frame(summary(mlm)$coefficients))
    }
               """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    df_stats["project"] = [val.replace("project", "") for val in df_stats.index.values]
    # print(df_stats)
    # df_stats = df_stats[df_stats["project"] != "(Intercept)"]  # removing the intercept line
    # df_stats = df_stats[["project", "Pr(>|z|)"]]
    return df_stats


def anova_gene_stats(dataframe, filename, name1, name2):
    """
    Perform an anova test on ``dataframe``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :param name1: (string) the name of a column in dataframe
    :param name2: (string) the
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    for colnames in dataframe.columns:
        if "size" in colnames:
            dataframe[colnames] = np.log10(dataframe[colnames].values)
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""

    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        mlm <- lm(%s ~ %s + project, data=data)
        png(paste(name, "_dignostics.png", sep=""), height=2160, width=870)
        par(mfrow=c(2, 2))
        plot(mlm)
        dev.off()
        mlm <- aov(%s ~ %s + project, data=data)
        return(as.data.frame(TukeyHSD(mlm, which="project")$project))
    }

                   """ % (name1, name2, name1, name2))
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    df_stats["project"] = df_stats.index.values
    df_stats = df_stats[["project", "diff", "lwr", "upr", "p adj"]]
    return df_stats
