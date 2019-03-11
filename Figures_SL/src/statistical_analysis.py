#!/usr/bin/env python3

# -*- coding: utf-8 -*-

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import figS3C_stat

pandas2ri.activate()


def anova_nt_stats(dataframe, filename):
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
    df_stats["project"] = [val.replace("project", "") for val in df_stats.index.values]
    df_stats = df_stats[df_stats["project"] != "(Intercept)"]  # removing the intercept line
    # df_stats = df_stats[["project", "Pr(>|t|)"]]
    return df_stats


def nb_glm_stats(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("MASS")
    require("DHARMa")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
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
    print(df_stats)
    df_stats = df_stats[df_stats["project"] != "(Intercept)"]  # removing the intercept line
    # df_stats = df_stats[["project", "Pr(>|z|)"]]
    return df_stats


def poisson_glm_stats(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        mlm <- glm(values ~ project, data=data, family="poisson")
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
    print(df_stats)
    df_stats = df_stats[df_stats["project"] != "(Intercept)"]  # removing the intercept line
    # df_stats = df_stats[["project", "Pr(>|z|)"]]
    return df_stats


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