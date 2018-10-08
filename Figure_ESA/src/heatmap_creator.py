#!/urs/bin/python3.5

# coding: utf-8

# imports
import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
import numpy as np
import figure_producer
import exon_control_handler
import union_dataset_function
import os
import random
import sys
sch = __import__('scipy.cluster.hierarchy')
nt_list = ["A", "C", "G", "T", "S", "W", "R", "Y"]
dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

# functions
def create_matrix(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for non_iupac characteristics of an exon set.

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
                final_value = float(median_obs - control_dic[target_columns[j]]) / control_dic[target_columns[j]] * 100
                projects_tab[j + 1].append([final_value])
                project_res.append(final_value)
            projects_tab[0].append(project_res)
    return projects_tab, project_names


def create_matrix_iupac(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
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
    new_targets = [target_columns[i].replace("iupac", "%s_nt" % nt)
                   for i in range(len(target_columns)) for nt in nt_list]
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
            for j in range(len(target_columns)):
                for nt in nt_list:
                    values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, target_columns[j], nt))
                    median_obs = np.median(values[~np.isnan(values)])
                    final_value = float(median_obs - control_dic[target_columns[j]][nt]) / \
                        control_dic[target_columns[j]][nt] * 100
                    project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def create_matrix_alldnt(cnx, id_projects, names, target_columns, control_dic, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for dnt characteristics of an exon set.

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
    new_targets = [target_columns[i].replace("dnt", "%s_dnt" % dnt)
                   for i in range(len(target_columns)) for dnt in remove_dnt(control_dic, target_columns[i])]
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
            for j in range(len(target_columns)):
                new_dnt_list = remove_dnt(control_dic, target_columns[j])
                for dnt in new_dnt_list:
                    values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, target_columns[j], dnt))
                    median_obs = np.median(values[~np.isnan(values)])
                    if control_dic[target_columns[j]][dnt] != 0:
                        final_value = float(median_obs - control_dic[target_columns[j]][dnt]) / \
                                      control_dic[target_columns[j]][dnt] * 100
                    else:
                        final_value = float("nan")
                    project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def create_matrix_nt(cnx, id_projects, names, target_column, control_dic, nt, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for a single nt of an exon set.

    :param cnx: (sqlite3 connect object) connexion to sed database
    :param id_projects: (list of ints) the list of splicing lore project id.
    :param names: (list of strings) the list of projects name (corresponding - in the same order - of the projects \
    in ``id_projects``) or if union is not none: list of sf_name.
    :param target_column: (string) characteristic of interest for a set of exons.
    :param control_dic: (dictionary of float) dictionary storing medians values for every characteristics of \
    teh control set of exons.
    :param nt: (string) a nucleotides.
    :param regulations: (list of strings) the strings can be "up" or "down" only for up or down-regulated exons.
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (lif of list of float) the medians value for a project (line) for the characteristic of \
    interest on one nucleotide (number of value in one line corresponding to a project).
    """
    new_targets = target_column.replace("iupac", "%s_nt" % nt)
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
            values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt))
            median_obs = np.median(values[~np.isnan(values)])
            final_value = float(median_obs - control_dic[target_column][nt]) / control_dic[target_column][nt] * 100
            project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def create_matrix_dnt(cnx, id_projects, names, target_column, control_dic, dnt, regulations, union=None):
    """
    Create a matrix of relative medians (toward control) for a single dnt of an exon set.

    :param cnx: (sqlite3 connect object) connexion to sed database
    :param id_projects: (list of ints) the list of splicing lore project id.
    :param names: (list of strings) the list of projects name (corresponding - in the same order - of the projects \
    in ``id_projects``) or if union is not none: list of sf_name.
    :param target_column: (string) characteristic of interest for a set of exons.
    :param control_dic: (dictionary of float) dictionary storing medians values for every characteristics of \
    teh control set of exons.
    :param dnt: (string) a di-nucleotides.
    :param regulations: (list of strings) the strings can be "up" or "down" only for up or down-regulated exons.
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (lif of list of float) the medians value for a project (line) for the characteristic of \
    interest on one nucleotide (number of value in one line corresponding to a project).
    """
    new_targets = target_column.replace("dnt", "%s_dnt" % dnt)
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            project_names.append("%s_%s" % (names[i], regulation))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
            values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, target_column, dnt))
            median_obs = np.median(values[~np.isnan(values)])
            if control_dic[target_column][dnt] != 0:
                final_value = float(median_obs - control_dic[target_column][dnt]) / control_dic[target_column][dnt] * 100
            else:
                final_value = float("nan")
            project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def simple_heatmap(data_array, labelsy, labelsx, output, name=""):
    """
    Create a heatmap link to a single dendrogram.

    :param data_array: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    :param labelsy: (list of strings) list of projects name
    :param labelsx: (list of strings) list of characteristics of interest
    :param output: (string)  path where the results will be created
    :param name: (string) partial name of the file
    """
    for i in range(len(data_array)):
        data_array[i][0] += random.random() / 1000000000
    data_array = np.array(data_array)
    labelsx = [labelsx]
    d = {labelsy[i]: {labelsx[j]: [i, j] for j in range(len(labelsx))} for i in range(len(labelsy))}
    # Initialize figure by creating upper dendrogram
    # Create Side Dendrogram
    figure = ff.create_dendrogram(data_array, orientation='right', labels=labelsy, linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'weighted'))
    for i in range(len(figure['data'])):
        figure['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure

    # Create Heatmap
    dendro_leaves = figure['layout']['yaxis']['ticktext']
    dendro_leaves2 = labelsx
    data_arrange = np.array([[None] * len(data_array[0])] * len(data_array))
    for i in range(len(dendro_leaves)):
        for j in range(len(dendro_leaves2)):
            vx = d[dendro_leaves[i]][dendro_leaves2[j]][0]
            vy = d[dendro_leaves[i]][dendro_leaves2[j]][1]
            data_arrange[i][j] = data_array[vx][vy]
    for i in range(len(data_arrange)):
        data_arrange[i][0] = round(data_arrange[i][0], 7)
    heatmap = [
        go.Heatmap(
            x=dendro_leaves2,
            y=dendro_leaves,
            z=data_arrange,
            colorbar={"x": -0.05},
            colorscale="Picnic",
            #colorscale = [[0.0, 'rgb(0, 114, 178)'], [0.25, 'rgb(86, 180, 233)'], [0.5, 'rgb(255, 255, 255)'], [0.75, 'rgb(240, 228, 66)'], [1.0, 'rgb(230, 159, 0)']],
            zmin=-20,
            zmax=20
        )
    ]
    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = figure['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    figure['data'].extend(heatmap)

    # Edit Layout
    figure['layout'].update({"title": "Clustering of '%s' for %s exons" % (labelsx[0], name),
                             "autosize": True, "height": 1080, "width": 1920,
                             'showlegend': False, 'hovermode': 'closest',
                             })
    # Edit xaxis
    figure['layout']['xaxis'].update({'domain': [0.805, 0.9],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks': "", })
    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .8],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': True,
                                      'ticks': "",
                                      "side": "right"})

    plotly.offline.plot(figure, filename='%s%s_%s.html' % (output, labelsx[0], name),
                        auto_open=False)


def heatmap_creator(data_array, labelsx, labelsy, output, name=""):
    """
    Create a heatmap link to two dendrograms.

    :param data_array: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    :param labelsy: (list of strings) list of projects name
    :param labelsx: (list of strings) list of characteristics of interest
    :param output: (string)  path where the results will be created
    :param name: (string) partial name of the file
    """
    for i in range(len(data_array)):
        data_array[i][0] += random.random() / 1000000000
    d = {labelsy[i]: {labelsx[j]: [i, j] for j in range(len(labelsx))} for i in range(len(labelsy))}
    # Initialize figure by creating upper dendrogram
    data_side = data_array.transpose()
    figure = ff.create_dendrogram(data_side, orientation='bottom', labels=labelsx, linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'weighted'))

    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array, orientation='right', labels=labelsy, linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'weighted'))
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    figure['data'].extend(dendro_side['data'])
    figure["layout"]["yaxis"] = dendro_side["layout"]["yaxis"]

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves2 = figure['layout']['xaxis']['ticktext']
    data_arrange = np.array([[None] * len(data_array[0])] * len(data_array))
    for i in range(len(dendro_leaves)):
        for j in range(len(dendro_leaves2)):
            vx = d[dendro_leaves[i]][dendro_leaves2[j]][0]
            vy = d[dendro_leaves[i]][dendro_leaves2[j]][1]
            data_arrange[i][j] = data_array[vx][vy]

    heatmap = [
        go.Heatmap(
            x=dendro_leaves2,
            y=dendro_leaves,
            z=data_arrange,
            colorbar={"x": -0.05},
            colorscale="Picnic",
            #colorscale=[[0.0, 'rgb(0, 114, 178)'], [0.25, 'rgb(86, 180, 233)'], [0.5, 'rgb(255, 255, 255)'],
            #            [0.75, 'rgb(240, 228, 66)'], [1.0, 'rgb(230, 159, 0)']],
            zmin=-20,
            zmax=20
        )
    ]

    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    figure['data'].extend(heatmap)

    # Edit Layout
    figure['layout'].update({"autosize": True, "height": 1080, "width": 1920,
                             'showlegend': False, 'hovermode': 'closest',
                             })
    # Edit xaxis
    figure['layout']['xaxis'].update({'domain': [0.15, 0.8],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks': ""})
    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0.11, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': True,
                                      'ticks': "",
                                      "side": "right"})
    # Edit yaxis2
    figure['layout'].update({'yaxis2': {'domain': [.825, 1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})
    plotly.offline.plot(figure, filename='%s%s.html' % (output, name),
                        auto_open=False)


def remove_dnt(control_dic, target_column):
    """
    Remove the dnt where control frequencies equals to 0.

    :param control_dic: (dictionary of float) dictionary storing medians values for every characteristics of \
    teh control set of exons.
    :param target_column: (string) column_of_interest
    :return: (list of string) list of dnt of interest
    """
    new_dnt_list = []
    for dnt in dnt_list:
        if control_dic[target_column][dnt] != 0:
            new_dnt_list.append(dnt)
    return new_dnt_list



def main():
    """
    Launch the main function.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = exon_control_handler.control_handler(cnx, exon_type)
    if len(sys.argv) < 2:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/heatmap/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = figure_producer.get_interest_project(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # Creating heatmap
            target_columns = ["gene_size", "nb_intron_gene", "median_intron_size", "upstream_exon_size", "exon_size",
                              "downstream_exon_size", "upstream_intron_size", "downstream_intron_size"]
            projects_tab, project_names = create_matrix(cnx, id_projects, name_projects, target_columns, ctrl_dic,
                                                        regulations)
            heatmap_creator(np.array(projects_tab[0]), target_columns, project_names, output,
                            "non_iupac_" + "_".join(regulations))
            # Creating dendrograms
            for i in range(1, len(projects_tab)):
                simple_heatmap(projects_tab[i], project_names, target_columns[i-1], output, "_".join(regulations))

            # handling iupac heatmap
            #  All iupac at once
            target_columns = ["iupac_exon", "iupac_gene", "iupac_upstream_intron_dist", "iupac_upstream_intron_proxi",
                              "iupac_downstream_intron_proxi", "iupac_downstream_intron_dist"]
            projects_tab, project_names, new_targets = create_matrix_iupac(cnx, id_projects, name_projects,
                                                                           target_columns, ctrl_dic, regulations)
            heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_iupac_" + "_".join(regulations))
            # Iupac by iupac
            for target_column in target_columns:
                projects_tab, project_names, new_targets = create_matrix_iupac(cnx, id_projects, name_projects,
                                                                               [target_column], ctrl_dic, regulations)
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                for nt in nt_list:
                    projects_tab, project_names, new_target = create_matrix_nt(cnx, id_projects, name_projects,
                                                                               target_column, ctrl_dic, nt, regulations)
                    simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

            # handling dnt heatmap
            #  All dnt at once
            target_columns = ["dnt_exon", "dnt_gene", "dnt_upstream_intron_dist", "dnt_upstream_intron_proxi",
                              "dnt_downstream_intron_proxi", "dnt_downstream_intron_dist"]
            projects_tab, project_names, new_targets = create_matrix_alldnt(cnx, id_projects, name_projects,
                                                                            target_columns, ctrl_dic, regulations)
            heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_dnt_" + "_".join(regulations))
            # dnt by dnt
            for target_column in target_columns:
                projects_tab, project_names, new_targets = create_matrix_alldnt(cnx, id_projects, name_projects,
                                                                                [target_column], ctrl_dic, regulations)
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                new_dnt_list = remove_dnt(ctrl_dic, target_column)
                for dnt in new_dnt_list:
                    projects_tab, project_names, new_target = create_matrix_dnt(cnx, id_projects, name_projects,
                                                                                target_column, ctrl_dic, dnt,
                                                                                regulations)
                    simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

    elif sys.argv[1] == "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/heatmap_union/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = union_dataset_function.get_splicing_factor_name(cnx)
        for regulations in [["up"], ["down"], ["up", "down"]]:
            # Creating heatmap
            target_columns = ["gene_size", "nb_intron_gene", "median_intron_size", "upstream_exon_size", "exon_size",
                              "downstream_exon_size", "upstream_intron_size", "downstream_intron_size"]
            projects_tab, project_names = create_matrix(cnx, None, name_projects, target_columns, ctrl_dic,
                                                        regulations, "union")
            heatmap_creator(np.array(projects_tab[0]), target_columns, project_names, output,
                            "non_iupac_" + "_".join(regulations))
            # Creating dendrograms
            for i in range(1, len(projects_tab)):
                simple_heatmap(projects_tab[i], project_names, target_columns[i - 1], output, "_".join(regulations))

            # handling iupac heatmap
            #  All iupac at once
            target_columns = ["iupac_exon", "iupac_gene", "iupac_upstream_intron_dist", "iupac_upstream_intron_proxi",
                              "iupac_downstream_intron_proxi", "iupac_downstream_intron_dist"]
            projects_tab, project_names, new_targets = create_matrix_iupac(cnx, None, name_projects,
                                                                           target_columns, ctrl_dic, regulations,
                                                                           "union")
            heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_iupac_" + "_".join(regulations))
            # Iupac by iupac
            for target_column in target_columns:
                projects_tab, project_names, new_targets = create_matrix_iupac(cnx, None, name_projects,
                                                                               [target_column], ctrl_dic, regulations,
                                                                               "union")
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                for nt in nt_list:
                    projects_tab, project_names, new_target = create_matrix_nt(cnx, None, name_projects,
                                                                               target_column, ctrl_dic, nt, regulations,
                                                                               "union")
                    simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

            # handling dnt heatmap
            #  All dnt at once
            target_columns = ["dnt_exon", "dnt_gene", "dnt_upstream_intron_dist", "dnt_upstream_intron_proxi",
                              "dnt_downstream_intron_proxi", "dnt_downstream_intron_dist"]
            projects_tab, project_names, new_targets = create_matrix_alldnt(cnx, None, name_projects,
                                                                            target_columns, ctrl_dic, regulations,
                                                                            "union")
            heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                            "all_dnt_" + "_".join(regulations))
            # dnt by dnt
            for target_column in target_columns:
                projects_tab, project_names, new_targets = create_matrix_alldnt(cnx, None, name_projects,
                                                                                [target_column], ctrl_dic, regulations,
                                                                                "union")
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output,
                                target_column + "_" + "_".join(regulations))
                new_dnt_list = remove_dnt(ctrl_dic, target_column)
                for dnt in new_dnt_list:
                    projects_tab, project_names, new_target = create_matrix_dnt(cnx, None, name_projects,
                                                                                target_column, ctrl_dic, dnt,
                                                                                regulations, "union")
                    simple_heatmap(projects_tab, project_names, new_target, output, "_".join(regulations))

if __name__ == "__main__":
    main()
