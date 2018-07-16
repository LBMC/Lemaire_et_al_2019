#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import heatmap_creator
import os
import figure_producer
import exon_control_handler
import sys
import numpy as np
import union_dataset_function

def main():
    """
    Launch the main function.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    if len(sys.argv) < 2:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/heatmap_proxi/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # handling iupac heatmap
            #  All iupac at once
            target_columns = ["iupac_upstream_intron_proxi",
                              "iupac_downstream_intron_proxi"]
            projects_tab, project_names, new_targets = heatmap_creator.create_matrix_iupac(cnx, id_projects, name_projects,
                                                                           target_columns, ctrl_dic, regulations)
            heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_iupac_" + "_".join(regulations))
            # Iupac by iupac
            for target_column in target_columns:
                projects_tab, project_names, new_targets = heatmap_creator.create_matrix_iupac(cnx, id_projects, name_projects,
                                                                               [target_column], ctrl_dic, regulations)
                heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                # for nt in heatmap_creator.nt_list:
                #     projects_tab, project_names, new_target = heatmap_creator.create_matrix_nt(cnx, id_projects, name_projects,
                #                                                                target_column, ctrl_dic, nt, regulations)
                #     heatmap_creator.simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

            # handling dnt heatmap
            #  All dnt at once
            target_columns = ["dnt_upstream_intron_proxi",
                              "dnt_downstream_intron_proxi"]
            projects_tab, project_names, new_targets = heatmap_creator.create_matrix_alldnt(cnx, id_projects, name_projects,
                                                                            target_columns, ctrl_dic, regulations)
            heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_dnt_" + "_".join(regulations))
            # dnt by dnt
            for target_column in target_columns:
                projects_tab, project_names, new_targets = heatmap_creator.create_matrix_alldnt(cnx, id_projects, name_projects,
                                                                                [target_column], ctrl_dic, regulations)
                heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                # for dnt in new_dnt_list:
                #     projects_tab, project_names, new_target = create_matrix_dnt(cnx, id_projects, name_projects,
                #                                                                 target_column, ctrl_dic, dnt,
                #                                                                 regulations)
                #     simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

    elif sys.argv[1] == "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/heatmap_proxi_union/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = union_dataset_function.get_splicing_factor_name(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # Creating heatmap
            target_columns = ["iupac_upstream_intron_proxi",
                              "iupac_downstream_intron_proxi"]
            # Iupac by iupac
            for target_column in target_columns:
                projects_tab, project_names, new_targets = heatmap_creator.create_matrix_iupac(cnx, None, name_projects,
                                                                               [target_column], ctrl_dic, regulations,
                                                                               "union")
                heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                # for nt in heatmap_creator.nt_list:
                #     projects_tab, project_names, new_target = heatmap_creator.create_matrix_nt(cnx, None, name_projects,
                #                                                                target_column, ctrl_dic, nt, regulations,
                #                                                                "union")
                #     heatmap_creator.simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

            # handling dnt heatmap
            #  All dnt at once
            target_columns = ["dnt_upstream_intron_proxi",
                              "dnt_downstream_intron_proxi"]
            projects_tab, project_names, new_targets = heatmap_creator.create_matrix_alldnt(cnx, None, name_projects,
                                                                            target_columns, ctrl_dic, regulations,
                                                                            "union")
            heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_dnt_" + "_".join(regulations))
            # dnt by dnt
            for target_column in target_columns:
                projects_tab, project_names, new_targets = heatmap_creator.create_matrix_alldnt(cnx, None, name_projects,
                                                                                [target_column], ctrl_dic, regulations,
                                                                                "union")
                heatmap_creator.heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                # new_dnt_list = heatmap_creator.remove_dnt(ctrl_dic, target_column)
                # for dnt in new_dnt_list:
                #     projects_tab, project_names, new_target = heatmap_creator.create_matrix_dnt(cnx, None, name_projects,
                #                                                                 target_column, ctrl_dic, dnt,
                #                                                                 regulations, "union")
                #     heatmap_creator.simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

if __name__ == "__main__":
    main()
