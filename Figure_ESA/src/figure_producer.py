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
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7, "K": 8, "M": 9}


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
    idp = [val[0] for val in res]
    name = [val[1] for val in res]
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
    for exon in exon_list:
        query = """SELECT %s
                   FROM sed
                   where gene_id = %s
                   AND exon_pos = %s """ % (target_column, exon[0], exon[1])
        cursor.execute(query)
        r = cursor.fetchone()[0]
        if r is not None:
            res.append(r)
    return res


def get_list_of_value_iupac(cnx, exon_list, target_column, nt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt: (string) a nucleotide
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
            res.append(float(r.split(";")[nt_dic[nt]]))
    return res


def get_values_for_many_projects(cnx, id_projects, target_column, regulation):
    """
    Return the value of ``target_column`` for each ``regulation`` exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string)) up or down
    :return: (list of list of float) each sublist of float corresponds to the values of ``target_column`` \
    for every regulated exon in a given project.
    """
    results = []
    for id_project in id_projects:
        exon_list = get_ase_events(cnx, id_project, regulation)
        results.append(get_list_of_value(cnx, exon_list, target_column))
    return results


def get_values_for_many_projects_iupac(cnx, id_projects, target_column, regulation, nt):
    """
    Return the frequency of the nucleotide ``nt`` of ``target_column`` for each ``regulation`` \
    exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param nt: (string) a nucleotide
    :return: (list of list of float) each sublist of float corresponds to the values of ``target_column`` \
    for every regulated exon in a given project.
    """
    results = []
    for id_project in id_projects:
        exon_list = get_ase_events(cnx, id_project, regulation)
        results.append(get_list_of_value_iupac(cnx, exon_list, target_column, nt))
    return results


def create_figure(cnx, id_projects, name_projects, target_column, regulation, output):
    """
    Create a figure for every column in sed database whose name does not contain "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param name_projects: (list of string) the list of project name in the same order that list id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    """
    result = get_values_for_many_projects(cnx, id_projects, target_column, regulation)
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
                         "box": {"visible": True}})
        else:
            data.append({"y": new_result[i], "type": "violin",
                         "name": new_name[i], "visible": 'legendonly', "marker": {"color": c[i]},
                         "box": {"visible": True}})
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


def create_figure_iupac(cnx, id_projects, name_projects, target_column, regulation, output, nt):
    """
    Create a figure for every column in sed database whose name contains "iupac".

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects: (list) list project id
    :param name_projects: (list of string) the list of project name in the same order that list id
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param output: (string) path where the result will be created
    :param nt: (string) the nt of interest
    """
    result = get_values_for_many_projects_iupac(cnx, id_projects, target_column, regulation, nt)
    target_column = target_column.replace("iupac", "%s_nt" % nt)
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
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/project_figures/"
    # If the output directory does not exist, then we create it !
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    regs = ["up", "down"]
    cnx = connexion(seddb)
    columns = get_column_of_interest(cnx)
    id_projects, name_projects = get_interest_project(cnx)
    for regulation in regs:
        print(regulation)
        for target_column in columns:
            print("   %s" % target_column)
            if "iupac" not in target_column:
                create_figure(cnx, id_projects, name_projects, target_column, regulation, output)
            else:
                for nt in nt_dic.keys():
                    create_figure_iupac(cnx, id_projects, name_projects, target_column, regulation, output, nt)

if __name__ == "__main__":
    main()
