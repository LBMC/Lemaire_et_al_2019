#!/usr/bin/python3.5

# coding : utf8

"""
Description:

    This script will display for each of up or down exon in a particular project every on of it's characteristics?
    This will use the sed database
"""


# import
import sqlite3
import plotly.graph_objs as go
import numpy as np
import plotly
import os
import sys
import union_dataset_function
import group_factor
import math
import control_exon_adapter
import function
import function_mfe
import exon_class
import exon_class_mfe
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
import rpy2.robjects.vectors as v
import pandas as pd
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7}
dnt_dic = {"AA": 0, "AC": 1, "AG": 2, "AT": 3, "CA": 4, "CC": 5,
           "CG": 6, "CT": 7, "GA": 8, "GC": 9, "GG": 10, "GT": 11,
           "TA": 12, "TC": 13, "TG": 14, "TT": 15}
log_columns = ["nb_intron_gene", "downstream_intron_size", "upstream_intron_size", "median_flanking_intron_size"]
exon_class.set_debug(0)
exon_class_mfe.set_debug(0)
size_bp_up_seq = 100
output_bp = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/bp_files/"


# Functions
def connexion(seddb):
    """
    Connexion to SED database.

    :param seddb: ((string) path to sed database
    :return:  (sqlite3 connection object) allow connexion to sed database
    """
    return sqlite3.connect(seddb)


def get_interest_project(cnx):
    """
    Get the id of every project defined in sed database (from splicing lore).

    :param cnx: (sqlite3 connection object) connexion to sed database
    :return: (list of int) list of id_project
    """
    cursor = cnx.cursor()
    query = "SELECT id, project_name FROM rnaseq_projects"
    cursor.execute(query)
    res = cursor.fetchall()
    idp = []
    name = []
    for val in res:
        if val[0] not in group_factor.bad_id_projects:
            idp.append(val[0])
            name.append(val[1])
    return idp, name


def get_ase_events(cnx, id_project, regulation):
    """
    Get every exon up or down regulated in a particular project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :param regulation: (string)) up or down
    :return: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    """
    if regulation == "up":
        regulation = ">= 0.1"
    else:
        regulation = "<= -0.1"
    cursor = cnx.cursor()
    query = """SELECT gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue_glm_cor <= 0.05""" % (id_project, regulation)
    cursor.execute(query)
    res = cursor.fetchall()
    if len(res) == 0:
            query = """SELECT gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue <= 0.05""" % (id_project, regulation)
            cursor.execute(query)
            res = cursor.fetchall()
    return res


def get_list_of_value(cnx, exon_list, target_column):
    """
    Get the individual values for ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["gene_size", "nb_intron_gene", "median_intron_size", "iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                res.append(r)
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                           FROM sed
                           where gene_id = %s
                           AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                if r is not None:
                    res.append(r)
                redundancy_gene_dic[exon[0]] = 1
    return res


def get_redundant_list_of_value(cnx, exon_list, target_column):
    """
    Get the individual values for ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    for exon in exon_list:
        query = """SELECT %s
                   FROM sed
                   where gene_id = %s
                   AND exon_pos = %s """ % (target_column, exon[0], exon[1])
        cursor.execute(query)
        r = cursor.fetchone()[0]
        if r is not None:
            res.append(r)
        else:
            res.append(None)
    return res


def get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt_dnt: (string) a nucleotide or di_nucleotide
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                if len(nt_dnt) == 1:
                    res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                else:
                    res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                if r is not None:
                    if len(nt_dnt) == 1:
                        res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                    else:
                        res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
                redundancy_gene_dic[exon[0]] = 1
    return res


def get_redundant_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt_dnt: (string) a nucleotide or di_nucleotide
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    for exon in exon_list:
        query = """SELECT %s
                   FROM sed
                   where gene_id = %s
                   AND exon_pos = %s """ % (target_column, exon[0], exon[1])
        cursor.execute(query)
        r = cursor.fetchone()[0]
        if r is not None:
            if len(nt_dnt) == 1:
                res.append(float(r.split(";")[nt_dic[nt_dnt]]))
            else:
                res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
        else:
            res.append(None)
    return res


def handle_nb_bp_recovering(cnx, exon_list, output, sf_name, regulation, target):
    """
    Recover nb bp for the exon list regulated by ``sf_name``
    :param cnx: (sqlite3 connect object) connection to fasterDB database
    :param exon_list: (list of 2 int) gene id + exon pos in gene
    :param output: (string) folder where the result will be created
    :param sf_name: (string) the name of the splicing factor studied
    :param regulation: (string) regulation up or down
    :param target: (string) the value we want to recover
    :return: (list of int) the list of good quality branch point
    """
    output_file = "%s%s_%s_%s_nt.py" % (output, sf_name, regulation, size_bp_up_seq)
    if not os.path.isfile(output_file):
        new_exon_list = [exon_class.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list]
        bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, sequence_list, ag_count_list, hbound_list = function.bp_ppt_calculator(new_exon_list, size_bp_up_seq)
        with open(output_file, "w") as bp_file:
            bp_file.write("bp_score=%s\n" % str(bp_score_list))
            bp_file.write("ppt_score=%s\n" % str(ppt_score_list))
            bp_file.write("nb_bp=%s\n" % str(nb_bp_list))
            bp_file.write("nb_good_bp=%s\n" % str(nb_good_bp_list))
            bp_file.write("bp_seq=%s\n" % str(sequence_list))
            bp_file.write("ag_count=%s\n" % str(ag_count_list))
            bp_file.write("hbound=%s\n" % str(hbound_list))
    else:
        sys.path.insert(0, output)
        mod = __import__(output_file.split("/")[-1].replace(".py", ""))
        nb_good_bp_list = mod.nb_good_bp
        hbound_list = mod.hbound
        ag_count_list = mod.ag_count
    if target == "nb_good_bp":
        return nb_good_bp_list
    elif target == "hbound":
        return hbound_list
    else:
        return ag_count_list


def handle_mfe_recovering(cnx, exon_list, output, sf_name, regulation, target):
    """
    Recover nb bp for the exon list regulated by ``sf_name``
    :param cnx: (sqlite3 connect object) connection to fasterDB database
    :param exon_list: (list of 2 int) gene id + exon pos in gene
    :param output: (string) folder where the result will be created
    :param sf_name: (string) the name of the splicing factor studied
    :param regulation: (string) regulation up or down
    :param target: (string) the value we want to recover
    :return: (list of int) the list of good quality branch point
    """
    output_file = "%s%s_%s_mfe.py" % (output, sf_name, regulation)
    if not os.path.isfile(output_file):
        new_exon_list = [exon_class_mfe.ExonClass(cnx, str(exon[0]), int(exon[0]), int(exon[1])) for exon in exon_list]
        mfe_3ss, mfe_5ss = function_mfe.mfe_calculator(new_exon_list)
        with open(output_file, "w") as bp_file:
            bp_file.write("mfe_3ss=%s\n" % str(mfe_3ss))
            bp_file.write("mfe_5ss=%s\n" % str(mfe_5ss))
    else:
        sys.path.insert(0, output)
        mod = __import__(output_file.split("/")[-1].replace(".py", ""))
        mfe_3ss = mod.mfe_3ss
        mfe_5ss = mod.mfe_5ss
    if target == "mfe_3ss":
        return mfe_3ss
    else:
        return mfe_5ss


def get_values_for_many_projects(cnx, cnx_fasterdb, id_projects_sf_names, target_column, regulation, output_bp_file, union):
    """
    Return the value of ``target_column`` for each ``regulation`` exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param cnx_fasterdb: (sqlite3 connection object) connexion to fasterdb database
    :param id_projects_sf_names: (list of str or  int) list project id if union is none. List of sf_name \
    else
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string)) up or down
    :param output_bp_file: (string) path where the bp files will be created
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (list of list of float) each sublist of float corresponds to the values of ``target_column`` \
    for every regulated exon in a given project.
    """
    results = []
    if not union:
        for id_project in id_projects_sf_names:
            exon_list = get_ase_events(cnx, id_project, regulation)
            if target_column == "median_flanking_intron_size":
                values1 = np.array(get_redundant_list_of_value(cnx, exon_list, "upstream_intron_size"), dtype=float)
                values2 = np.array(get_redundant_list_of_value(cnx, exon_list, "downstream_intron_size"),dtype=float)
                values = np.array([np.nanmedian([values1[i], values2[i]]) for i in range(len(values1))])
                results.append(values)
            elif target_column in ["nb_good_bp", "hbound", "ag_count"]:
                results.append(handle_nb_bp_recovering(cnx_fasterdb, exon_list, output_bp_file, str(id_project), regulation,
                                                       target_column))
            elif "mfe" in target_column:
                results.append(handle_mfe_recovering(cnx_fasterdb, exon_list, output_bp_file, str(id_project), regulation, target_column))
            else:
                results.append(get_list_of_value(cnx, exon_list, target_column))

    else:
        for sf_name in id_projects_sf_names:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
            if target_column == "median_flanking_intron_size":
                values1 = np.array(get_redundant_list_of_value(cnx, exon_list, "upstream_intron_size"), dtype=float)
                values2 = np.array(get_redundant_list_of_value(cnx, exon_list, "downstream_intron_size"),dtype=float)
                values = np.array([np.nanmedian([values1[i], values2[i]]) for i in range(len(values1))])
                results.append(values)
            elif target_column == "min_flanking_intron_size":
                values1 = np.array(get_redundant_list_of_value(cnx, exon_list, "upstream_intron_size"), dtype=float)
                values2 = np.array(get_redundant_list_of_value(cnx, exon_list, "downstream_intron_size"),
                                   dtype=float)
                values = np.array([np.nanmin([values1[i], values2[i]]) for i in range(len(values1))])
                results.append(values)
            elif target_column in ["nb_good_bp", "hbound", "ag_count"]:
                results.append(handle_nb_bp_recovering(cnx_fasterdb, exon_list, output_bp_file, sf_name, regulation, target_column))
            elif "mfe" in target_column:
                results.append(handle_mfe_recovering(cnx_fasterdb, exon_list, output_bp_file, sf_name, regulation, target_column))
            else:
                results.append(get_list_of_value(cnx, exon_list, target_column))
    return results


def get_values_for_many_projects_iupac_dnt(cnx, id_projects_sf_names, target_column, regulation, nt_dnt, union):
    """
    Return the frequency of the nucleotide ``nt`` of ``target_column`` for each ``regulation`` \
    exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects_sf_names: (list of str or  int) list project id if union is none. List of sf_name \
    else
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param nt_dnt: (string) a nucleotide or a di-nucleotide
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (list of list of float) each sublist of float corresponds to the values of ``target_column`` \
    for every regulated exon in a given project.
    """

    results = []
    if not union:
        for id_project in id_projects_sf_names:
            exon_list = get_ase_events(cnx, id_project, regulation)
            results.append(get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt))

    else:
        for sf_name in id_projects_sf_names:
            exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
            results.append(get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt))
    return results


def mann_withney_test_r(list_values1, list_values2):
    """
    Perform a mann withney wilcoxon test on ``list_values1`` and ``list_values2``.

    :param list_values1: (list of float)  list of float
    :param list_values2: (list of float)  list of float
    :return: (float) the pvalue of the mann withney test done one `list_values1`` and ``list_values2``.
    """
    wicox = robj.r("""

    function(x, y){
        test = wilcox.test(x,y, alternative='two.sided', correct=F)
        return(test$p.value)
    }

                   """)
    pval = float(wicox(v.FloatVector(list_values1), v.FloatVector(list_values2))[0])
    return pval


def create_statistical_report(list_values, list_name, ctrl_full, filename):
    """
    Create a statistical report.

    :param list_values: (list of list of floats) the list of value that we want to compare to a control list
    :param list_name: (list of string) the name of each sublist of float in ``list_values``
    :param ctrl_full: (list of float) the control list of values
    :param output: (string) the file where the report will be created
    """
    cur_ctrl = np.array(ctrl_full, dtype=float)
    cur_ctrl = list(cur_ctrl[~np.isnan(cur_ctrl)])
    dic_res = {"Factor": [], "P-value": []}
    for i in range(len(list_values)):
        dic_res["Factor"].append(list_name[i])
        cur_list = np.array(list_values[i], dtype=float)
        cur_list = list(cur_list[~np.isnan(cur_list)])
        # print("          Factor : %s, mean = %s" % (list_name[i], np.nanmean(list_values[i])))

        dic_res["P-value"].append(mann_withney_test_r(cur_list, cur_ctrl))
    df = pd.DataFrame(dic_res)
    rstats = robj.packages.importr('stats')
    pcor = rstats.p_adjust(v.FloatVector(dic_res["P-value"]), method="BH")
    df["P-adjusted_BH"] = pcor
    df.to_csv(filename.replace(".html", "_stat.txt"), sep="\t", index=False)


def create_figure(cnx, cnx_fasterdb, id_projects, name_projects, target_column, regulation, output, ctrl_dic, output_bp_file, ctrl_full, union=None):
    """
    Create a figure for every column in sed database whose name does not contain "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param cnx_fasterdb: (sqlite3 connection object) connexion to fasterdb database
    :param id_projects: (list) list project id
    :param name_projects: (list of string) the list of project name (or sf name) in the same order that list id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    :param ctrl_dic: (string) control dictionary
    :param output_bp_file: (string) path where the bp files will be created
    :param ctrl_full: (string) full control dictionary
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    """
    if not union:
        result = get_values_for_many_projects(cnx, cnx_fasterdb, id_projects, target_column, regulation, output_bp_file, union)
    else:
        result = get_values_for_many_projects(cnx, cnx_fasterdb, name_projects, target_column, regulation, output_bp_file, union)
    d = {name_projects[i]: result[i] for i in range(len(result))}
    e = sorted(d.items(), key=lambda x: np.median(x[1]), reverse=True)
    new_name = [x[0] for x in e]
    new_result_tmp = [x[1] for x in e]
    log = False
    if target_column in log_columns:
        log = True
    if not log:
        new_result = [x[1] for x in e]
    else:
        new_result = [list(map(math.log10, x[1])) for x in e]

    if not log:
        new_result = [(np.array(val) - ctrl_dic[target_column]) / ctrl_dic[target_column] * 100 for val in new_result]
        title = '%s of %s exons for every projects' % (target_column, regulation)
    else:
        title = 'log10 %s of %s exons for every projects' % (target_column, regulation)
        new_result = [(np.array(val) - math.log10(ctrl_dic[target_column])) / math.log10(ctrl_dic[target_column]) * 100 for val in new_result]

    cb = ['hsl(' + str(h) + ',50%' + ',60%)' for h in np.linspace(0, 360, len(new_result))]
    cv = ['hsl(' + str(h) + ',50%' + ',80%)' for h in np.linspace(0, 360, len(new_result))]
    data = []
    for i in range(len(new_result)):
        data.append({"y": new_result[i], "type": "violin",
                     "name": new_name[i], "visible": True, "fillcolor": cv[i], "opacity": 1, "line": {"color": "black"},
                     "box": {"visible": True, "fillcolor" : cb[i]}, "meanline": {"visible": False}})
    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(200, 200, 200)',
            gridwidth=1,
            zerolinecolor='rgb(200, 0, 0)',
            zerolinewidth=2,
            title="relative %s" % target_column
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
    filename = "%s%s_%s_exons_figure.html" % (output, target_column, regulation)
    plotly.offline.plot(fig, filename=filename,
                        auto_open=False, validate=False)
    create_statistical_report(new_result_tmp, new_name, ctrl_full[target_column], filename)


def create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, nt_dnt, ctrl_dic, ctrl_full, union=None):
    """
    Create a figure for every column in sed database whose name contains "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param name_projects: (list of string) the list of project name (or sf_name if union is not none) \
     in the same order that list id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    :param nt_dnt: (string) the nt of interest or the di-nucleotide of interest
    :param ctrl_dic: (string) control dictionary
    :param ctrl_full: (string) full control dictionary
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    """
    if not union:
        result = get_values_for_many_projects_iupac_dnt(cnx, id_projects, target_column, regulation, nt_dnt, union)
    else:
        result = get_values_for_many_projects_iupac_dnt(cnx, name_projects, target_column, regulation, nt_dnt, union)
    d = {name_projects[i]: result[i] for i in range(len(result))}
    e = sorted(d.items(), key=lambda x: np.median(x[1]), reverse=True)
    new_name = [x[0] for x in e]
    new_result_tmp = [x[1] for x in e]
    #mean_val = np.nanmean(ctrl_full[target_column][nt_dnt])
    #new_result = [(np.array(val) - mean_val) / mean_val * 100 for val in new_result]
    new_result = [(np.array(val) - ctrl_dic[target_column][nt_dnt]) / ctrl_dic[target_column][nt_dnt] * 100 for val in new_result_tmp]
    cb = ['hsl(' + str(h) + ',50%' + ',60%)' for h in np.linspace(0, 360, len(new_result))]
    cv = ['hsl(' + str(h) + ',50%' + ',80%)' for h in np.linspace(0, 360, len(new_result))]
    target_column_new = target_column.replace("iupac", "%s_nt" % nt_dnt)
    target_column_new = target_column_new.replace("dnt", "%s_dnt" % nt_dnt)
    data = []
    for i in range(len(new_result)):
        data.append({"y": new_result[i], "type": "violin",
                     "name": new_name[i], "visible": True, "fillcolor": cv[i], "opacity": 1, "line": {"color": "black"},
                     "box": {"visible": True, "fillcolor" : cb[i]}, "meanline": {"visible": False}})
    layout = go.Layout(
        title='%s of %s exons for every projects' % (target_column_new, regulation),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(200, 200, 200)',
            gridwidth=1,
            zerolinecolor='rgb(200, 0, 0)',
            zerolinewidth=2,
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
    filename = "%s%s_%s_exons_figure.html" % (output, target_column_new, regulation)
    plotly.offline.plot(fig, filename=filename,
                        auto_open=False, validate=False)
    # print("    Stat %s" % (target_column_new))
    create_statistical_report(new_result_tmp, new_name, ctrl_full[target_column][nt_dnt], filename)


def main():
    """
    Launch the creation of figures.
    """
    exon_type = "CCE"
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    fasterdb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/fasterDB_lite.db"
    regs = ["up", "down"]
    cnx = connexion(seddb)
    cnx_fasterdb = connexion(fasterdb)
    #columns = ["iupac_exon", "exon_size", "upstream_intron_size", "downstream_intron_size", "gene_size", "median_flanking_intron_size", "force_donor", "force_acceptor", "iupac_upstream_intron_adjacent1", "nb_intron_gene", "nb_good_bp_%s" % size_bp_up_seq, "hbound", "ag_count", "mfe_3ss", "mfe_5ss", "iupac_upstream_intron_ppt_area"]
    #columns = ["force_donor", "force_acceptor", "upstream_intron_size", "downstream_intron_size","iupac_exon", "iupac_upstream_intron_proxi", "iupac_downstream_intron_proxi", "min_flanking_intron_size"]
    columns = ["iupac_exon"]
    ctrl_dic, ctrl_full = control_exon_adapter.control_handler(cnx, exon_type, size_bp_up_seq)
    if len(sys.argv) < 2:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/project_figures_new/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = get_interest_project(cnx)
        for regulation in regs:
            print(regulation)
            for target_column in columns:
                if target_column == "nb_good_bp":
                    if not os.path.isdir(output_bp):
                        os.mkdir(output_bp)
                print("   %s" % target_column)
                if "iupac" in target_column:
                    for nt in nt_dic.keys():
                        create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, nt, ctrl_dic, ctrl_full)
                elif "dnt" in target_column:
                    for dnt in dnt_dic.keys():
                        create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, dnt, ctrl_dic, ctrl_full)
                else:
                    create_figure(cnx, cnx_fasterdb, id_projects, name_projects, target_column, regulation, output, ctrl_dic, output_bp, ctrl_full)
    elif sys.argv[1] == "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/project_figures_union_new/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = group_factor.get_wanted_sf_name(cnx)
        for regulation in regs:
            print(regulation)
            for target_column in columns:
                if target_column == "nb_good_bp":
                    if not os.path.isdir(output_bp):
                        os.mkdir(output_bp)
                print("   %s" % target_column)
                if "iupac" in target_column:
                    for nt in nt_dic.keys():
                        create_figure_iupac_dnt(cnx, None, name_projects, target_column, regulation, output, nt, ctrl_dic, ctrl_full,
                                                "union")
                elif "dnt" in target_column:
                    for dnt in dnt_dic.keys():
                        create_figure_iupac_dnt(cnx, None, name_projects, target_column, regulation, output, dnt, ctrl_dic, ctrl_full,
                                                "union")
                else:
                    create_figure(cnx, cnx_fasterdb, None, name_projects, target_column, regulation, output, ctrl_dic, output_bp, ctrl_full, "union")
    else:
        print("wrong arg !")
    cnx.close()
    cnx_fasterdb.close()

if __name__ == "__main__":
    main()
