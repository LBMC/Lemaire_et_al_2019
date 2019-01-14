#!/usr/bin/python3

"""
Description:
    The goal of this script is to create an heatmap with each line corresponds to an exon set found in a file \
    contained in a given folder
"""

import heatmap_maker
import numpy as np
import figure_producer
import os
import control_exon_adapter
import argparse
import pandas as pd
import plotly
import plotly.graph_objs as go


def extract_exon(filename):
    """
    Extract exon from a file ``filename``.

    :param filename: (string) a file containing exons
    :return: (list of list of 2 int) list of exon identified by their gene id and exon_name
    """
    exon_list = []
    with open(filename, "r") as my_file:
        count = 0
        for line in my_file:
            count += 1
            line = line.replace("\n", "").split("_")
            try:
                gene_id = int(line[0])
                exon_pos = int(line[1])
                exon_list.append([gene_id, exon_pos])
            except ValueError:
                print("Line %s, does not contain an exon" % count)
    return exon_list


def create_matrix(cnx, list_files, names_files, target_columns, control_dic):
    """
    Create a matrix of relative medians (toward control) for iupac characteristics of an exon set.

    :param cnx: (sqlite3 connect object) connexion to sed database
    :param list_files: (list of string)  the list of file of interest
    :param names_files: (list of strings) the list of file name
    :param target_columns: (list of strings) list of interest characteristics for a set of exons.
    :param control_dic: (dictionary of float) dictionary storing medians values for every characteristics of \
    teh control set of exons.
    :return: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    """
    new_targets = heatmap_maker.create_columns_names(target_columns)
    project_names = []
    projects_tab = []
    for i in range(len(names_files)):
        project_res = []
        project_names.append("%s" % (names_files[i]))
        exon_list = extract_exon(list_files[i])
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
                else:

                    values = np.array(figure_producer.get_list_of_value(cnx, exon_list, new_targets[j]))
                median_obs = np.median(values[~np.isnan(values)])
                final_value = float(median_obs - control_dic[new_targets[j]]) / \
                    control_dic[new_targets[j]] * 100
            project_res.append(final_value)
        projects_tab.append(project_res)
    return projects_tab, project_names, new_targets


def heatmap_gc_sorted(data_array, labelsx, labelsy, output, contrast, name=""):
    """
    Create a GC sorted heatmap, sorted on gc content if they are presents in labels X

    :param data_array: (lif of list of float) the medians value for a project (line) for every characteristic of \
    interests (number of value in one line corresponding to a project).
    :param labelsy: (list of strings) list of projects name
    :param labelsx: (list of strings) list of characteristics of interest
    :param output: (string)  path where the results will be created
    :param contrast: (int) the value of the contrast
    :param name: (string) partial name of the file
    """
    df = pd.DataFrame(data_array, index=labelsy)
    index_g = []
    index_c = []
    index_s = []
    for i in range(len(labelsx)):
        if "S" in labelsx[i]:
            index_s.append(i)
        if "G" in labelsx[i]:
            index_g.append(i)
        if "C" in labelsx[i]:
            index_c.append(i)
    if len(index_s) > 0:
        final_df = df.sort_values(index_s[0], ascending=True)
    elif len(index_g) > 0 and len(index_c) > 0:
        df["m"] = (df[index_g[0]] +  df[index_c[0]]) / 2
        df = df.sort_values("m", ascending=True)
        final_df = df.drop("m", axis=1)
    else:
        final_df = df.sort_values(0, ascending=True)
    if final_df is not None:
        heatmap = [go.Heatmap(z=final_df.values,
                              x=labelsx,
                              y=list(final_df.index),
                              colorbar={"x": -0.05},
                              colorscale="Picnic",
                              zmin=-contrast,
                              zmax=contrast)]
        layout = go.Layout(yaxis = dict(
                                          autorange=True,
                                          ticks="-",
                                          side="right"),
                            xaxis = dict(ticks=""),
            margin=dict(t=200, r=500, b=200, l=200)
        )
        figure = go.Figure(data=heatmap, layout=layout)
        plotly.offline.plot(figure, filename='%s%s_sorted.html' % (output, name),
                        auto_open=False)

def main(columns, name, contrast, input_fodler):
    """
    Launch the main function.
    """
    tmp = os.listdir(input_fodler)
    list_file = ["%s/%s" % (input_fodler, tmp[i]) for i in range(len(tmp))]
    list_name = [tmp[i].replace(".txt", "") for i in range(len(tmp))]
    target_columns = columns.split(",")
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = figure_producer.connexion(seddb)
    ctrl_dic = control_exon_adapter.control_handler(cnx, exon_type)
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/heatmap_from_file_exon_list/"
    # If the output directory does not exist, then we create it !
    if not os.path.isdir(output):
        os.mkdir(output)
    for regulations in [["down"]]:
        # Creating heatmap

        projects_tab, project_names, new_targets = create_matrix(cnx, list_file, list_name,
                                                                 target_columns, ctrl_dic)
        # if len(new_targets) > 1:
        #     heatmap_maker.heatmap_creator(np.array(projects_tab), new_targets, project_names, output, contrast, name)
        # else:
        #     heatmap_maker.simple_heatmap(np.array(projects_tab), new_targets, project_names, output, contrast, name)
        heatmap_gc_sorted(np.array(projects_tab), new_targets, project_names, output, contrast, name)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Create a heatmap with  the wanted parameters

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("required argument")


    required_args.add_argument('--columns', dest='columns',
                               help="The columns (in sed database) you want to analyze, they must be coma separated",
                               required=True)

    required_args.add_argument('--name', dest='name',
                               help="""name of heatmap""",
                               required=True)
    required_args.add_argument('--input', dest='input',
                               help="""A folder containing files contaiing an exon list""",
                               required=True)

    parser.add_argument("--nt", dest="nt", help="the list ow wanted nucleoitdes coma separted",
                        default="A,T,G,C")
    parser.add_argument("--contrast", dest="contrast", help="value for the color-scale of the heatmap",
                        default=50)
    args = parser.parse_args()  # parsing arguments

    global nt_list
    nt_list = args.nt.split(",")
    heatmap_maker.make_global(nt_list)
    for nt in nt_list:
        if nt not in ["A", "C", "G", "T", "S", "W", "K", "M"]:
            parser.error("Wrong value given in the parameter nt")
    try:
        args.contrast = int(args.contrast)
    except ValueError:
        parser.error("contrast argument must be an integer !")
    main(args.columns, args.name, args.contrast, args.input)  # executing the program


if __name__ == "__main__":
    launcher()
