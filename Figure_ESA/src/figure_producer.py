#!/usr/bin/python3.5

# coding : utf8

"""
Description:

    This script will display for each of up or down exon in a particular project every on of it's characteristics?
    This will use the sed database
"""


# import
import sqlite3


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
    return cursor.fetchall()


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
        res.append(cursor.fetchone()[0])
    return res

