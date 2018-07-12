#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# imports
import numpy as np
import figure_producer
import control_exon_adapter
import union_dataset_function
import os
import sys
import heatmap_creator
sch = __import__('scipy.cluster.hierarchy')


# Functions
def get_control_information(exon_type, control_file):
    """
    Test whether we have all the information we need for ``exon_type`` exons

    :param exon_type: (string) the type of control exon we want to use
    :param control_file: (string) the file were we want to write the control data results
    :return: None if the control file does'nt exist, or a dictionary that contains \
    all the information about the CCE exon
    """
    control_folder = os.path.dirname(control_file)
    if not os.path.isdir(control_folder):
        os.mkdir(control_folder)
        return None
    if not os.path.isfile(control_file):
        return None
    sys.path.insert(0, control_folder)
    mod = __import__("control")
    try:
        dic = eval("mod.%s_dic" % exon_type)
    except AttributeError:
        return None
    if "rel_exon_size" not in dic.keys():
        print("Adapting")
    return dic


def create_matrix(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for relative_size characteristic of an exon set. \

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
    relative_size_dic = {"rel_exon_intron_up": ["exon_size", "upstream_intron_size"],
                         "rel_exon_intron_down" : ["exon_size", "downstream_intron_size"],
                         "rel_exon_introns": ["exon_size", "upstream_intron_size", "downstream_intron_size"]}
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
                if len(relative_size_dic[target_columns[j]]) == 2:
                    target1 = relative_size_dic[target_columns[j]][0]
                    target2 = relative_size_dic[target_columns[j]][1]
                    values1 = np.array(figure_producer.get_list_of_value(cnx, exon_list, target1))
                    values2 = np.array(figure_producer.get_list_of_value(cnx, exon_list, target2))
                    relative_values = np.array(values1, dtype=float) / np.array(values2, dtype=float)

                elif len(relative_size_dic[target_columns[j]]) == 3:
                    target1 = relative_size_dic[target_columns[j]][0]
                    target2 = relative_size_dic[target_columns[j]][1]
                    target3 = relative_size_dic[target_columns[j]][2]
                    values1 = np.array(figure_producer.get_list_of_value(cnx, exon_list, target1))
                    values2 = np.array(figure_producer.get_list_of_value(cnx, exon_list, target2))
                    values3 = np.array(figure_producer.get_list_of_value(cnx, exon_list, target3))
                    relative_values = np.array(values1, dtype=float) / ((np.array(values2, dtype=float) + np.array(values3, dtype=float)) / 2)
                else:
                    print("something went wrong ! ")
                    exit(1)
                tmp = relative_values[~np.isnan(relative_values)]
                median_obs = np.median(tmp[~np.isinf(tmp)])
                final_value = float(median_obs - control_dic[target_columns[j]]) / control_dic[target_columns[j]] * 100
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
    ctrl_dic = control_exon_adapter.control_handler(cnx, exon_type)
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/relative_size_heatmap/"
    target_columns = ["rel_exon_intron_up", "rel_exon_intron_down", "rel_exon_introns"]
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
                                            "relative_size_" + "_".join(regulations))
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
