#!/usr/bin/env python3

# -*- coding: utf-8 -*-

from rpy2.robjects import r, pandas2ri

pandas2ri.activate()


def glm_poisson_stats(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    dataframe = dataframe[dataframe["project"] != "CCE_exons"]
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require("multcomp")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        data$project <- as.factor(as.vector(data$project))
        mlm <- glm(values ~ project, data=data, family="poisson")
        simulationOutput <- simulateResiduals(fittedModel = mlm, n = 250)
        png(paste(name, "_dignostics.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        s <- summary(glht(mlm, mcp(project = "Tukey")))
        tab <- as.data.frame(cbind(names(s$test$coefficients), s$test$coefficients, s$test$sigma,  s$test$tstat, s$test$pvalues))
        colnames(tab) <- c("Comparison", "Estimate", "Std.Error", "z.value", "Pr(>|z|)")
        return(tab)
    }
               """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    return df_stats


def glm_poisson_stats_spliceosome(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    dataframe = dataframe[(dataframe["project"] != "CCE_exons") & (dataframe["project"] != "CCE")]
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require("multcomp")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        data$project <- as.factor(as.vector(data$project))
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
    df_stats["project"] = df_stats.index
    return df_stats


def glm_nb_stats(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    dataframe = dataframe[dataframe["project"] != "CCE_exons"]
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require("multcomp")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        data$project <- as.factor(as.vector(data$project))
        mlm <- glm.nb(values ~ project, data=data)
        simulationOutput <- simulateResiduals(fittedModel = mlm, n = 250)
        png(paste(name, "_dignostics.png", sep=""), height=1080, width=1920)
        par(mfrow=c(2, 2))
        plot(simulationOutput)
        dev.off()
        s <- summary(glht(mlm, mcp(project = "Tukey")))
        tab <- as.data.frame(cbind(names(s$test$coefficients), s$test$coefficients, s$test$sigma,  s$test$tstat, s$test$pvalues))
        colnames(tab) <- c("Comparison", "Estimate", "Std.Error", "z.value", "Pr(>|z|)")
        return(tab)
    }
               """)
    name = filename.split(".")[0]
    df_stats = pandas2ri.ri2py(stat_s(r_df, name))
    return df_stats


def glm_nb_stats_spliceosome(dataframe, filename):
    """
    Perform a glm nb test on ``list_values1`` and ``list_values2``.

    :param dataframe: (pandas DataFrame) a dataframe
    :param filename: (string)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    dataframe = dataframe[(dataframe["project"] != "CCE_exons") & (dataframe["project"] != "CCE")]
    r_df = pandas2ri.py2ri(dataframe)
    stat_s = r("""
    require("DHARMa")
    require("multcomp")
    function(data, name){
        factors <- levels(as.factor(as.vector(data$project)))
        data$project <- as.factor(as.vector(data$project))
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
    df_stats["project"] = df_stats.index
    return df_stats