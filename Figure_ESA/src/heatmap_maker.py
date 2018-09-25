#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is to create very flexible heatmap with every parameters of interest
"""


import pandas as pd
import numpy as np
import figure_producer
import union_dataset_function
import heatmap_creator
import exon_control_handler
import argparse
import os
import sys

nt_list = ["A", "T", "G", "C"]

def create_columns_names(target_columns):
    """
    Create columns name list.
    :param target_columns: (list of strings) list of columns names in sed

    :return: (list of string) list of columns names in the heatmap (if iupac is mentionned each nucleotides must \
    be specified.
    """
    new_targets = []
    for target in target_columns:
        if "iupac" in target:
            new_targets += [target.replace("iupac", "%s_nt" % nt) for nt in nt_list]
        else:
            new_targets.append(target)
    return new_targets


def create_matrix(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for iupac characteristics of an exon set.

    :param cnx: (sqlite3 connect object) connexion to sed database
    :param id_projects: (list of ints) the list of splicing lore project id.
    :param names: (list of strings) the list of projects name (corresponding - in the same order - of the projects \
    in ``id_projects``) or if union is not none: list of sf_name.
    :param target_columns: (list of strings) list of interest characteristics for a set of exons.
    :param control_dic: (dictionary of float) dictionary storing medians values for every characteristics of \
    teh control set of exons.
    :param regulations: (list of strings) the strings can be "up" or "down" only for up or down-regulated exons.
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    """
    new_targets = create_columns_names(target_columns)
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_events_4_a_sl(cnx, names[i], regulation)
            for j in range(len(new_targets)):
                if "_nt_" in new_targets[j]:
                    nt = new_targets[j].split("_")[0]
                    name_col = new_targets[j].replace("%s_nt" % nt, "iupac")
                    values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, name_col, nt))
                    median_obs = np.median(values[~np.isnan(values)])
                    final_value = float(median_obs - control_dic[name_col][nt]) / \
                        control_dic[name_col][nt] * 100
                else:
                    values = np.array(figure_producer.get_list_of_value(cnx, exon_list, new_targets[j]))
                    median_obs = np.median(values[~np.isnan(values)])
                    final_value = float(median_obs - control_dic[new_targets[j]]) / \
                        control_dic[new_targets[j]] * 100
                project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def main(union, columns, name):
    """
    Launch the main function.
    """
    target_columns = columns.split(",")
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    if union != "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/new_heatmap/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects_full, name_projects_full = figure_producer.get_interest_project(cnx)
        id_projects = []
        name_projects = []
        for i in range(len(name_projects_full)):
            if "SRRM4" not in name_projects_full[i] and "DNMT3B" not in name_projects_full[i]:
                name_projects.append(name_projects_full[i])
                id_projects.append(id_projects_full[i])
        for regulations in [["down"]]:
            # Creating heatmap

            projects_tab, project_names, new_targets = create_matrix(cnx, id_projects, name_projects,
                                                                           target_columns, ctrl_dic, regulations)
            heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            name)


    else:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/new_heatmap_union/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects_full = union_dataset_function.get_splicing_factor_name(cnx)
        name_projects = [name for name in name_projects_full if name != "DNMT3B" and name != "SRRM4"]
        name_projects = ["SNRNP70", "SNRPC", "SF1", "SF3A3", "SF3B1", "SF3B4", "U2AF1", "U2AF2"]
        for regulations in [["down"]]:
            # Creating heatmap

            projects_tab, project_names, new_targets = create_matrix(cnx, None, name_projects,
                                                                           target_columns, ctrl_dic, regulations,
                                                                           "union")
            heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output, name)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Create a heatmap with  the wanted parameters

    """,
                                     usage='%(prog)s --union union --columns columns')
    # Arguments for the parser
    required_args = parser.add_argument_group("required argument")

    parser.add_argument('--union', dest='union',
                               help="""union if you want create an heatmap on union dataset, nothing else""",
                              default="")

    required_args.add_argument('--columns', dest='columns',
                               help="""The columns (in sed database) you want to analyze, they must be coma separated""",
                              required=True)

    required_args.add_argument('--name', dest='name',
                               help="""name of heatmap""",
                              required=True)
    args = parser.parse_args()  # parsing arguments


    main(args.union, args.columns, args.name)  # executing the program


if __name__ == "__main__":
    launcher()