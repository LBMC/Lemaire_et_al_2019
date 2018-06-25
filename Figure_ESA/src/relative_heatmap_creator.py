#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# imports
import numpy as np
import figure_producer
import exon_control_handler
import union_dataset_function
import os
import sys
import heatmap_creator
sch = __import__('scipy.cluster.hierarchy')


# Functions

def create_matrix(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for non_iupac characteristics of an exon set. \
    The characteristics "such as force_donor and force acceptor" will be corrected using the control dictionary. \
    The relative frequencies of donor and acceptor will not be corrected because they are already corrected.

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
    interests (number of value in one line corresponding to a project). If `relative` not in target column, \
    the relative median is corrected with the control value in the control dictionary.
    """
    project_names = []
    projects_tab = [[] for i in range(len(target_columns) + 1)]
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
            for j in range(len(target_columns)):
                values = np.array(figure_producer.get_list_of_value(cnx, exon_list, target_columns[j]))
                median_obs = np.median(values[~np.isnan(values)])
                if "relative" in target_columns[j]:
                    final_value = median_obs
                else:
                    final_value = float(median_obs - control_dic[target_columns[j]]) / control_dic[
                        target_columns[j]] * 100
                projects_tab[j + 1].append([final_value])
                project_res.append(final_value)
            projects_tab[0].append(project_res)
    return projects_tab, project_names


def main():
    """
    Launch the main function.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/relative_heatmap/"
    target_columns = ["force_donor", "force_acceptor", "relative_donor_upstream", "relative_donor_downstream",
                      "relative_acceptor_upstream", "relative_acceptor_downstream"]
    if len(sys.argv) < 2:
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # Creating heatmap
            projects_tab, project_names = create_matrix(cnx, id_projects, name_projects, target_columns, ctrl_dic,
                                                        regulations)
            heatmap_creator.heatmap_creator(np.array(projects_tab[0]), target_columns, project_names, output,
                                            "forces_" + "_".join(regulations))

    elif sys.argv[1] == "union":
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = union_dataset_function.get_splicing_factor_name(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # Creating heatmap
            projects_tab, project_names = create_matrix(cnx, None, name_projects, target_columns, ctrl_dic,
                                                        regulations, "union")
            heatmap_creator.heatmap_creator(np.array(projects_tab[0]), target_columns, project_names, output,
                                            "union_forces_" + "_".join(regulations))

if __name__ == "__main__":
    main()
