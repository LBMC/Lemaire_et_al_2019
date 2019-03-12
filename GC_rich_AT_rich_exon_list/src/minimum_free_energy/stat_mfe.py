#!/usr/bin/env python3

# -*- coding: utf-8 -*-

from rpy2.robjects import r, pandas2ri

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
        png(paste(name, "_distrib_qqplots.png", sep=""), height=2160, width=1920)
        par(mfrow = c(3, 3))
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

        mlm <- aov(values ~ project, data=data)
        png(paste(name, "_dignostics.png", sep=""), height=2160, width=870)
        par(mfrow=c(4, 1))
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


def anova_nt_stats_spliceosome(dataframe, filename):
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
        png(paste(name, "_distrib_qqplots.png", sep=""), height=2160, width=1920)
        par(mfrow = c(3, 3))
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
    df_stats["project"] = df_stats.index.values
    # df_stats = df_stats[["project", "diff", "lwr", "upr", "p adj"]]
    return df_stats