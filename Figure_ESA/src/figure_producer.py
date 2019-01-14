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
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7, "K": 8, "M": 9}
dnt_dic = {"AA": 0, "AC": 1, "AG": 2, "AT": 3, "CA": 4, "CC": 5,
           "CG": 6, "CT": 7, "GA": 8, "GC": 9, "GG": 10, "GT": 11,
           "TA": 12, "TC": 13, "TG": 14, "TT": 15}


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


def get_values_for_many_projects(cnx, id_projects_sf_names, target_column, regulation, union):
    """
    Return the value of ``target_column`` for each ``regulation`` exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects_sf_names: (list of str or  int) list project id if union is none. List of sf_name \
    else
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string)) up or down
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


def create_figure(cnx, id_projects, name_projects, target_column, regulation, output, union=None):
    """
    Create a figure for every column in sed database whose name does not contain "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param name_projects: (list of string) the list of project name (or sf name) in the same order that list id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    """
    if not union:
        result = get_values_for_many_projects(cnx, id_projects, target_column, regulation, union)
    else:
        result = get_values_for_many_projects(cnx, name_projects, target_column, regulation, union)
    d = {name_projects[i]: result[i] for i in range(len(result))}
    e = sorted(d.items(), key=lambda x: np.median(x[1]), reverse=False)
    new_name = [x[0] for x in e]
    new_result = [x[1] for x in e]
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(new_result))]
    data = []
    for i in range(len(new_result)):
        if i in list(range(5)) + list(range(len(new_result)-5, len(new_result))):
            data.append({"y": new_result[i], "type": "violin",
                         "name": new_name[i], "visible": True, "marker": {"color": c[i]},
                         "box": {"visible": True}, "meanline": {"visible": True}})
        else:
            data.append({"y": new_result[i], "type": "violin",
                         "name": new_name[i], "visible": 'legendonly', "marker": {"color": c[i]},
                         "box": {"visible": True}, "meanline": {"visible": True}})
    layout = go.Layout(
        title='%s of %s exons for every projects' % (target_column, regulation),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_exons_figure.html" % (output, target_column, regulation),
                        auto_open=False, validate=False)


def create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, nt_dnt, union=None):
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
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    """
    if not union:
        result = get_values_for_many_projects_iupac_dnt(cnx, id_projects, target_column, regulation, nt_dnt, union)
    else:
        result = get_values_for_many_projects_iupac_dnt(cnx, name_projects, target_column, regulation, nt_dnt, union)
    target_column = target_column.replace("iupac", "%s_nt" % nt_dnt)
    target_column = target_column.replace("dnt", "%s_dnt" % nt_dnt)
    d = {name_projects[i]: result[i] for i in range(len(result))}
    e = sorted(d.items(), key=lambda x: np.median(x[1]), reverse=False)
    new_name = [x[0] for x in e]
    new_result = [x[1] for x in e]
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(new_result))]
    data = []
    for i in range(len(new_result)):
        if i in list(range(5)) + list(range(len(new_result)-5, len(new_result))):
            data.append({"y": new_result[i], "type": "violin",
                         "name": new_name[i], "visible": True, "marker": {"color": c[i]},
                         "box": {"visible": True}, "meanline": {"visible": True}})
        else:
            data.append({"y": new_result[i], "type": "violin",
                         "name": new_name[i], "visible": 'legendonly', "marker": {"color": c[i]},
                         "box": {"visible": True}, "meanline": {"visible": True}})
    layout = go.Layout(
        title='%s of %s exons for every projects' % (target_column, regulation),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            autotick=True,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
        ),
        margin=dict(
            l=40,
            r=30,
            b=150,
            t=100,
        ),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
    )

    fig = {"data": data, "layout": layout}
    plotly.offline.plot(fig, filename="%s%s_%s_exons_figure.html" % (output, target_column, regulation),
                        auto_open=False, validate=False)


def get_column_of_interest(cnx):
    """
    Get the columns name of interest in sed database.

    :param cnx: (sqlite3 object) connexion to sed database
    :return: (list of string) the name of the column of interest in sed table of the sed database.
    """
    cursor = cnx.cursor()
    cursor.execute("PRAGMA table_info(sed);")
    res = cursor.fetchall()
    columns = [a[1] for a in res]
    return columns[4:]


def main():
    """
    Launch the creation of figures.
    """
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    regs = ["up", "down"]
    cnx = connexion(seddb)
    columns = get_column_of_interest(cnx)
    if len(sys.argv) < 2:
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/project_figures/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        id_projects, name_projects = get_interest_project(cnx)
        for regulation in regs:
            print(regulation)
            for target_column in columns:
                print("   %s" % target_column)
                if "iupac" in target_column:
                    for nt in nt_dic.keys():
                        create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, nt)
                elif "dnt" in target_column:
                    for dnt in dnt_dic.keys():
                        create_figure_iupac_dnt(cnx, id_projects, name_projects, target_column, regulation, output, dnt)
                else:
                    create_figure(cnx, id_projects, name_projects, target_column, regulation, output)
    elif sys.argv[1] == "union":
        output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/project_figures_union/"
        # If the output directory does not exist, then we create it !
        if not os.path.isdir(output):
            os.mkdir(output)
        name_projects = union_dataset_function.get_splicing_factor_name(cnx)
        for regulation in regs:
            print(regulation)
            for target_column in columns:
                print("   %s" % target_column)
                if "iupac" in target_column:
                    for nt in nt_dic.keys():
                        create_figure_iupac_dnt(cnx, None, name_projects, target_column, regulation, output, nt,
                                                "union")
                elif "dnt" in target_column:
                    for dnt in dnt_dic.keys():
                        create_figure_iupac_dnt(cnx, None, name_projects, target_column, regulation, output, dnt,
                                                "union")
                else:
                    create_figure(cnx, None, name_projects, target_column, regulation, output, "union")
    else:
        print("wrong arg !")

if __name__ == "__main__":
    main()
