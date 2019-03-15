#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is to create very flexible heatmap with every parameters of interest
"""


import numpy as np
import figure_producer
import union_dataset_function
import control_exon_adapter
import argparse
import os
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
import random
import group_factor
import math
sch = __import__('scipy.cluster.hierarchy')


def redundant_ag_at_and_u1_u2(cnx, regulation):
    """
    Create the list of redundant exons between the AT and GC rich list of exons and \
    between the U1 and U2 list of exons

    :param cnx: (sqlite3 connect object) allow connection to sed database
    :param regulation: (string) the regulation we want for the common exons
    :return: (list of list of 2 int) list of exons identified by their gene id and their exons position
    """
    exon_at = []
    for sf_name in group_factor.at_rich_down:
        exon_at += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_at = union_dataset_function.washing_events_all(exon_at)
    exon_gc = []
    for sf_name in group_factor.gc_rich_down:
        exon_gc += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_gc = union_dataset_function.washing_events_all(exon_gc)
    global redundant_gc_at
    redundant_gc_at = [exon for exon in exon_at if exon in exon_gc]
    print("redundant exon GC and AT rich : %s" % len(redundant_gc_at))

    exon_u1 = []
    for sf_name in group_factor.u1_factors:
        exon_u1 += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_u1 = union_dataset_function.washing_events_all(exon_u1)
    exon_u2 = []
    for sf_name in group_factor.u2_factors:
        exon_u2 += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_u2 = union_dataset_function.washing_events_all(exon_u2)
    global redundant_u1_u2
    redundant_u1_u2 = [exon for exon in exon_u1 if exon in exon_u2]
    print("redundant exon U1 and U2 rich : %s" % len(redundant_u1_u2))


def difference(cnx, exon_list1, name_exon_list, regulation, sf_type):
    """
    Remove the redundant exons within the ``exon_list1``  i.e the exons regulated by both u1 and u2 \
    if ``name_exon_list`` corresponds to a component a the spliceosome or the exons regulated by \
    both AT and GC factors if ``name_exon_list`` corresponds to a splicing factors name.
    :param cnx: (sqlite3 connection object) connection to sed database
    :param exon_list1: (list of list of 2 int) a list of exons (each exon are identified by its gene id and its \
    position within the gene) regulated in ``name_exon_list`` project
    :param name_exon_list: (string) a splicing factor or a project's name.
    :param regulation: (string) the regulation of each exon within the exon list
    :param sf_type: (string) the type of splicing factor we want to display in the final figures
    :return: (list of list of 2 int) the ``exon_list1`` without the the exons regulated by both u1 and u2 \
    if ``name_exon_list`` corresponds to a component a the spliceosome or the exons regulated by \
    both AT and GC factors if ``name_exon_list`` corresponds to a splicing factors name.
    """
    if sf_type in ["GC_rich", "AT_rich", "spliceosome"]:
        print("----------------- removing redundant exon of %s " % name_exon_list)
        if "_" in name_exon_list:
            name_exon_list = name_exon_list.split("_")[0]
        if name_exon_list in group_factor.u1_factors + group_factor.u2_factors:
            exon_list_pur = [exon for exon in exon_list1 if exon not in redundant_u1_u2]
            print("-------> %s exon before, %s after" % (len(exon_list1), len(exon_list_pur)))
            return exon_list_pur
        elif name_exon_list in group_factor.gc_rich_down + group_factor.at_rich_down:
            exon_list_pur = [exon for exon in exon_list1 if exon not in redundant_gc_at]
            print("-------> %s exon before, %s after" % (len(exon_list1), len(exon_list_pur)))
            return exon_list_pur
        else:
            return exon_list1
    else:
        return exon_list1


def create_columns_names(target_columns):
    """
    Create columns name list.
    :param target_columns: (list of strings) list of columns names in sed
    :return: (list of string) list of columns names in the heatmap (if iupac is mentionned each nucleotides must \
    be specified.
    """
    new_targets = []
    for target in target_columns:
        if "iupac_upstream_intron_adjacent1" == target:
            new_targets += [target.replace("iupac", "Y_nt")]
        elif "iupac" in target:
            new_targets += [target.replace("iupac", "%s_nt" % nt) for nt in nt_list]
        else:
            new_targets.append(target)
    return new_targets


def create_matrix(cnx, id_projects, names, target_columns, control_dic, regulations, union=None, sf_type=None):
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
    :param sf_type: (string) the type of splicing factor we want to display in the final figures
    :return: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    """
    new_targets = create_columns_names(target_columns)
    project_names = []
    projects_tab = []
    for i in range(len(names)):
        for regulation in regulations:
            project_res = []
            if len(regulations) > 1:
                project_names.append("%s_%s" % (names[i], regulation))
            else:
                project_names.append("%s" % (names[i]))
            if not union:
                exon_list = figure_producer.get_ase_events(cnx, id_projects[i], regulation)
                print("Splicing factor : %s, project %s  - exons %s" % (names[i], id_projects[i], len(exon_list)))
                exon_list = difference(cnx, exon_list, names[i], regulation, sf_type)
            else:
                exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, names[i], regulation)
                print("Splicing factor : %s - exons %s" % (names[i], len(exon_list)))
                exon_list = difference(cnx, exon_list, names[i], regulation, sf_type)
            for j in range(len(new_targets)):
                if "_nt_" in new_targets[j]:
                    nt = new_targets[j].split("_")[0]
                    name_col = new_targets[j].replace("%s_nt" % nt, "iupac")
                    if "mean_intron" in new_targets[j]:
                        values1 = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, "iupac_upstream_intron", nt))
                        values2 = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, "iupac_downstream_intron", nt))
                        values = np.array([np.nanmedian([values1[i], values2[i]]) for i in range(len(values1))])
                    else:
                        values = np.array(figure_producer.get_list_of_value_iupac_dnt(cnx, exon_list, name_col, nt))

                    median_obs = np.median(values[~np.isnan(values)])
                    final_value = float(median_obs - control_dic[name_col][nt]) / \
                        control_dic[name_col][nt] * 100
                else:
                    if new_targets[j] == "median_flanking_intron_size":
                        values1 = np.array(
                            figure_producer.get_redundant_list_of_value(cnx, exon_list, "upstream_intron_size"),
                            dtype=float)
                        values2 = np.array(
                            figure_producer.get_redundant_list_of_value(cnx, exon_list, "downstream_intron_size"),
                            dtype=float)
                        values = np.array([np.nanmedian([values1[i], values2[i]]) for i in range(len(values1))])
                    elif new_targets[j] in ["nb_good_bp", "hbound", "ag_count"]:
                        if union:
                            values = np.array(figure_producer.handle_nb_bp_recovering(cnx_fasterdb, exon_list, output_bp, names[i], regulation, new_targets[j]))
                        else:
                            values = np.array(figure_producer.handle_nb_bp_recovering(cnx_fasterdb, exon_list, output_bp, str(id_projects[i]), regulation, new_targets[j]))
                    elif "mfe" in new_targets[j]:
                        if union:
                            values = np.array(figure_producer.handle_mfe_recovering(cnx_fasterdb, exon_list, output_bp, names[i], regulation, new_targets[j]))
                        else:
                            values = np.array(figure_producer.handle_mfe_recovering(cnx_fasterdb, exon_list, output_bp, str(id_projects[i]), regulation, new_targets[j]))
                    else:
                        values = np.array(figure_producer.get_list_of_value(cnx, exon_list, new_targets[j]))
                    if new_targets[j] in figure_producer.log_columns:
                        median_obs = np.median(values[~np.isnan(values)])
                        final_value = float(math.log10(median_obs) - math.log10(control_dic[new_targets[j]])) / math.log10(control_dic[
                            new_targets[j]]) * 100
                    else:
                        median_obs = np.median(values[~np.isnan(values)])
                        final_value = float(median_obs - control_dic[new_targets[j]]) / control_dic[new_targets[j]] * 100
                project_res.append(final_value)
            projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def project_t2_remove(cnx, sf_list):
    """
    Get the id of the projects we want to remove
    :param cnx: (sqlite3 connect object) connection to sed database
    :param sf_list: (string) list of splicing factor
    :return: (list of int) the id project to remove
    """
    cursor = cnx.cursor()
    list_projects = []
    query = "SELECT id FROM rnaseq_projects WHERE sf_name = '%s'"
    for sf in sf_list:
        cursor.execute(query % sf)
        res = cursor.fetchall()
        for project in res:
            list_projects.append(project[0])
    return list_projects


def removing_projects(full_projects_list, name_projects_list, project_2_remove):
    """
    Return ``full_projects_list`` without every project present in ``project_2_remove``

    :param full_projects_list: (list of int) list of every project in splicing lore
    :param name_projects_list:  (list of string) the nam of every project in splicing lore
    :param project_2_remove: (list of int) list of every project we need to remove in splicing lore
    :return: (2 lists):
        * (list of int) the list of project we want to analye
        * (list of string) the list of projects name we want to  analyse
    """
    wanted_project_id = []
    wanted_project_name = []
    for i in range(len(full_projects_list)):
        if full_projects_list[i] not in project_2_remove:
            wanted_project_id.append(full_projects_list[i])
            wanted_project_name.append(name_projects_list[i])
    return wanted_project_id, wanted_project_name


def simple_heatmap(data_array, labelsx, labelsy, output, name=""):
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
    d = {labelsy[i]: {labelsx[j]: [i, j] for j in range(len(labelsx))} for i in range(len(labelsy))}
    # Initialize figure by creating upper dendrogram
    # Create Side Dendrogram
    figure = ff.create_dendrogram(data_array, orientation='right', labels=labelsy, linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'average'))
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
            zmin=-50,
            zmax=50
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

    print('%s%s_%s.html' % (output, labelsx[0], name))
    plotly.offline.plot(figure, filename='%s%s.html' % (output, name),
                        auto_open=False)


def heatmap_creator(data_array, labelsx, labelsy, output, name=""):
    """
    Create a heatmap.

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
    figure = ff.create_dendrogram(data_side, orientation='bottom', labels=labelsx,
                                   linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'average'))

    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array, orientation='right', labels=labelsy,
                                  linkagefun=lambda x: sch.cluster.hierarchy.linkage(x, 'average'))
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
            # colorscale=[[0.0, 'rgb(0, 114, 178)'], [0.25, 'rgb(86, 180, 233)'], [0.5, 'rgb(255, 255, 255)'],
            #            [0.75, 'rgb(240, 228, 66)'], [1.0, 'rgb(230, 159, 0)']],
            zmin=-50,
            zmax=50
        )
    ]
    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']
    #
    # # Add Heatmap Data to Figure
    figure['data'].extend(heatmap)

    # figure = go.Figure(data=heatmap)

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
                                      'showticklabels': True,
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


def heatmap_gc_sorted(data_array, labelsx, labelsy, output, name=""):
    """
    Create a GC sorted heatmap, sorted on gc content if they are presents in labels X

    :param data_array: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    :param labelsy: (list of strings) list of projects name
    :param labelsx: (list of strings) list of characteristics of interest
    :param output: (string)  path where the results will be created
    :param name: (string) partial name of the file
    """
    df = pd.DataFrame(data_array, index=labelsy)
    index_g = []
    index_c = []
    index_s = []
    for i in range(len(labelsx)):
        if "G" in labelsx[i]:
            index_g.append(i)
        if "C" in labelsx[i]:
            index_c.append(i)
        if "S" in labelsx[i]:
            index_s.append(i)
    if len(index_s) > 0:
        final_df = df.sort_values(index_s[0], ascending=False)
    elif len(index_g) > 0 and len(index_c) > 0:
        df["m"] = (df[index_g[0]] +  df[index_c[0]]) / 2
        df = df.sort_values("m", ascending=False)
        final_df = df.drop("m", axis=1)
    else:
        final_df = df.sort_values(0, ascending=True)
    if final_df is not None:
        heatmap = [go.Heatmap(z=final_df.values,
                   x=labelsx,
                   y=list(final_df.index),
                   colorbar={"x": -0.05},
                   colorscale="Picnic",
                   zmin=-50,
                   zmax=50)]
        layout = go.Layout(yaxis = dict(
                                          ticks="",
                                          side="right"),
                            xaxis = dict(ticks=""))
        figure = go.Figure(data=heatmap, layout=layout)
        plotly.offline.plot(figure, filename='%s%s_sorted.html' % (output, name),
                        auto_open=False)



def main(union, columns, name, sf_type, regulation):
    """
    Launch the main function.
    """
    target_columns = columns.split(",")
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = control_exon_adapter.control_handler(cnx, exon_type)
    global output_bp
    output_bp = figure_producer.output_bp
    fasterdb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/fasterDB_lite.db"
    global cnx_fasterdb
    cnx_fasterdb = figure_producer.connexion(fasterdb)
    print("test")
    if union != "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/new_heatmap/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = group_factor.get_id_and_name_project_wanted(cnx, sf_type)
        print(id_projects)
        print(name_projects)
        for regulations in [regulation]:
            # Creating heatmap
            if sf_type is not None:
                redundant_ag_at_and_u1_u2(cnx, regulations[0])
            projects_tab, project_names, new_targets = create_matrix(cnx, id_projects, name_projects,
                                                                     target_columns, ctrl_dic, regulations, None,
                                                                     sf_type)
            if len(new_targets) > 1:
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output, name)
            heatmap_gc_sorted(np.array(projects_tab), new_targets, project_names, output, name=name)

    else:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/new_heatmap_union/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = group_factor.get_wanted_sf_name(sf_type)
        for regulations in [regulation]:
            # Creating heatmap
            if sf_type is not None:
                redundant_ag_at_and_u1_u2(cnx, regulations[0])
            projects_tab, project_names, new_targets = create_matrix(cnx, None, name_projects,
                                                                     target_columns, ctrl_dic, regulations,
                                                                     "union", sf_type)
            print(project_names)
            if len(new_targets) > 1:
                heatmap_creator(np.array(projects_tab), new_targets, project_names, output, name)
            heatmap_gc_sorted(np.array(projects_tab), new_targets, project_names, output, name)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Create a heatmap with  the wanted parameters

    """,
                                     usage='%(prog)s --union union --columns columns --name [--sf_type --nt]')
    # Arguments for the parser
    required_args = parser.add_argument_group("required argument")

    parser.add_argument('--union', dest='union',
                        help="""union if you want create an heatmap on union dataset, nothing else""",
                        default="")

    required_args.add_argument('--columns', dest='columns',
                               help="The columns (in sed database) you want to analyze, they must be coma separated",
                               required=True)

    required_args.add_argument('--name', dest='name',
                               help="""name of heatmap""",
                               required=True)

    parser.add_argument("--sf_type", dest="sf_type", help="the type of sf we want to display on "
                                                          "the heatmap (GC_rich, AT_rich, spliceosome)", default=None)

    parser.add_argument("--nt", dest="nt", help="the list ow wanted nucleoitdes coma separted",
                        default="A,T,G,C")
    parser.add_argument("--reg", dest="reg", help="the regulation_wanted",
                        default="down")
    args = parser.parse_args()  # parsing arguments

    global nt_list
    nt_list = args.nt.split(",")
    for nt in nt_list:
        if nt not in ["A", "C", "G", "T", "S", "W", "K", "M", "R", "Y"]:
            parser.error("Wrong value given in the parameter nt")
    if args.sf_type == "None":
        args.sf_type = None
    if args.sf_type not in ["GC_rich", "AT_rich", "spliceosome", "CF", None]:
        parser.error("Wrong value given for argument sf_type")
    args.reg = args.reg.split(",")
    for my_reg in args.reg:
        if my_reg not in ["up", "down"]:
            parser.error("Wrong regulation ! ")
    main(args.union, args.columns, args.name, args.sf_type, args.reg)  # executing the program


if __name__ == "__main__":
    launcher()
