#!/usr/bin/env python3

"""
Description:
    displays for every given list  of exons their violin plots for a \
    wanted feature
"""


import lazyparser as lp
import os
import sys
import numpy as np
import plotly.graph_objs as go
import plotly
import math
base1 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, base1)
import group_factor
base2 = os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
base2 += "/Figures_SL/src"
sys.path.insert(0, base2)
import figure_producer as fp
import control_exon_adapter as cea


def read_bed(bed_file):
    """
    Get every line of a bed file.

    :param bed_file: (str) a bed file
    :return: (list of list of data) every line in the bed file
    """
    bed_list = []
    with open(bed_file, "r") as bedin:
        for line in bedin:
            line = line.replace("\n", "")
            line = line.split("\t")
            line[0] = int(line[0])
            line[1] = int(line[1])
            bed_list.append(line)
    return bed_list


def get_values_for_many_projects(cnx, cnx_fasterdb, exon_lists,
                                 target_column, output_bp_file):
    """
    Return the value of ``target_column`` for each ` exons for \
    list of of exons ``exon_lists``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param cnx_fasterdb: (sqlite3 connection object) connexion to fasterdb \
    database
    :param exon_lists: (list of str) list of files containing exons
    :param target_column: (string) the column for which we want to get \
    information on exons.
    :param output_bp_file: (string) path where the tmp files will be created
    :return: (list of list of float) each sublist of float corresponds to \
    the values of ``target_column``for exons in ``exon_lists``.
    """
    results = []
    for exon_file in exon_lists:
        exon_list = read_bed(exon_file)
        if target_column == "median_flanking_intron_size":
            values1 = np.array(
                fp.get_redundant_list_of_value(cnx, exon_list,
                                               "upstream_intron_size"),
                dtype=float)
            values2 = np.array(
                fp.get_redundant_list_of_value(cnx, exon_list,
                                               "downstream_intron_size"),
                dtype=float)
            values = np.array([np.nanmedian([values1[i], values2[i]])
                               for i in range(len(values1))])
            results.append(values)
        elif target_column == "min_flanking_intron_size":
            values1 = np.array(
                fp.get_redundant_list_of_value(cnx, exon_list,
                                               "upstream_intron_size"),
                dtype=float)
            values2 = np.array(
                fp.get_redundant_list_of_value(cnx, exon_list,
                                               "downstream_intron_size"),
                               dtype=float)
            values = np.array([np.nanmin([values1[i], values2[i]])
                               for i in range(len(values1))])
            results.append(values)
        elif target_column in ["nb_good_bp", "hbound", "ag_count"]:
            results.append(
                fp.handle_nb_bp_recovering(cnx_fasterdb, exon_list,
                                           output_bp_file, exon_list,
                                           "", target_column))
        elif "mfe" in target_column:
            results.append(
                fp.handle_mfe_recovering(cnx_fasterdb, exon_list,
                                         output_bp_file, exon_list,
                                         "", target_column))
        else:
            results.append(fp.get_list_of_value(cnx, exon_list, target_column))
    return results


def get_values_for_many_projects_iupac_dnt(cnx, exons_files, target_column,
                                           nt_dnt):
    """
    Return the frequency of the nucleotide ``nt`` of ``target_column`` \
    for each ``regulation`` exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exons_files: (list of str) list of exons file
    :param target_column: (string) the column for which we want to get \
    information on exons.
    :param nt_dnt: (string) a nucleotide or a di-nucleotide
    :return: (list of list of float) each sublist of float corresponds to \
    the values of ``target_column`` for every regulated exon in exon_list.
    """
    results = []
    for exon_file in exons_files:
        exon_list = read_bed(exon_file)
        results.append(fp.get_list_of_value_iupac_dnt(cnx, exon_list,
                                                      target_column,
                                                      nt_dnt))
    return results


def create_figure_iupac_dnt(cnx, cnx_fasterdb, list_file,  name_file,
                            target_column, output, ctrl_full, exon_type,
                            tmp_folder, nt_dnt=None):
    """
    Create a figure for every column in sed database.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param cnx_fasterdb: (sqlite3 connection object) connexion to fasterdb \
    database
    :param list_file: (list of string) list of  exon files.
    :param name_file: (list of string) the list of name of each exon files \
    given in ``list_file``
    :param target_column: (string) the column for which we want to get\
     information on exons.
    :param output: (string) path where the result will be created
    :param nt_dnt: (string) the nt of interest or the di-nucleotide of interest
    :param ctrl_full: (string) full control dictionary
    :param tmp_folder: (str) location where temporary files are stored.
    :param exon_type: (str) the type of control exon
    """
    log_columns = ["upstream_intron_size", "downstream_intron_size",
                   "min_flanking_intron_size", "exon_size"]
    if nt_dnt is not None:
        target_column_new = target_column.replace("iupac", "%s_nt" % nt_dnt)\
            .replace("dnt", "%s_dnt" % nt_dnt)
        result = get_values_for_many_projects_iupac_dnt(cnx, list_file,
                                                        target_column, nt_dnt)
    else:
        target_column_new = target_column
        result = get_values_for_many_projects(cnx, cnx_fasterdb, list_file,
                                              target_column, tmp_folder)
        if target_column in log_columns:
            result = [list(map(math.log10, res)) for res in result]
            target_column_new = "log10_" + target_column_new

    filename = "%s/%s_figure.html" % (output, target_column_new)
    # Adding control values to the ones of interest
    name_file += [exon_type]
    if nt_dnt is not None:
        my_ctrl = np.array(ctrl_full[target_column][nt_dnt], dtype=float)
    else:
        my_ctrl = np.array(ctrl_full[target_column], dtype=float)
    my_ctrl = list(my_ctrl[~np.isnan(my_ctrl)])
    if target_column not in log_columns:
        result += [my_ctrl]
    else:
        result += [list(map(math.log10, my_ctrl))]
    color_dic = group_factor.color_dic
    color_bright = group_factor.color_dic_bright
    default_colors = [color_dic["GC_pure"], color_dic["AT_pure"], "#F3A431",
                      "#A000A0", "#5c0001", color_dic["CCE"]]
    default_bright = [color_bright["GC_pure"], color_bright["AT_pure"],
                      "#FFDB65", "#d055d0", "#b05544", color_bright["CCE"]]
    if len(default_colors) != len(result):
        default_colors = default_colors[0:len(result) - 1] + \
            [default_colors[-1]]
        default_bright = default_bright[0:len(result) - 1] + \
            [default_bright[-1]]
    data = []
    for i in range(len(result)):
        data.append({"y": result[i], "type": "violin",
                     "name": name_file[i], "visible": True,
                     "fillcolor": default_bright[i],
                     "opacity": 1, "line": {"color": "black"},
                     "box": {"visible": True, "fillcolor": default_colors[i]},
                     "meanline": {"visible": False}})
    layout = go.Layout(
        title='%s' % target_column_new,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            # autotick=True,
            gridcolor='rgb(200, 200, 200)',
            gridwidth=1,
            title="relative %s" % target_column_new
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename=filename,
                        auto_open=False, validate=False)


@lp.parse(list_file="file", seddb="file", fasterdb="file",
          exon_type=["CCE"], output="dir")
def figure_creator(list_file, name_file, seddb, fasterdb, output, targets,
                   exon_type="CCE"):
    """
    Create violin plot figures representing a target parameter \
    for every list of exons given.

    :param list_file: (List(vtype=str)) List of exons file
    :param name_file: (List(vtype=str)) List of the name for the exon files.
    :param seddb: (str) File corresponding to sed database
    :param fasterdb: (str) File corresponding to fasterdb database
    :param output: (str) the folder where the result will be created
    :param targets: (List(vtype=str)) list of targets
    :param exon_type: (str) the type of control exons of interest
    """
    tmp_folder = output + "/tmp"
    if not os.path.isdir(tmp_folder):
        os.mkdir(tmp_folder)
    cnx = fp.connexion(seddb)
    cnx_fasterdb = fp.connexion(fasterdb)
    size = 100
    nt_list = ["S", "R"]
    ctrl_dic, ctrl_full = cea.control_handler(cnx, exon_type, size)
    for target in targets:
        print("Working on %s" % target)
        if "iupac" in target:
            for nt in nt_list:
                create_figure_iupac_dnt(cnx, cnx_fasterdb, list_file,
                                        name_file, target, output, ctrl_full,
                                        exon_type, tmp_folder, nt_dnt=nt)
        else:
            create_figure_iupac_dnt(cnx, cnx_fasterdb, list_file,
                                    name_file, target, output, ctrl_full,
                                    exon_type, tmp_folder, nt_dnt=None)
    cnx.close()
    cnx_fasterdb.close()


if __name__ == "__main__":
    figure_creator()
