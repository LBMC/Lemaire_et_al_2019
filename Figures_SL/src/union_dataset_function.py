#!/usr/bin/env python3

# coding: utf-8

import group_factor
import numpy as np


def get_gene_name(cnx, gene_id):
    """
    Give the gene name thanks to a sedDB gene id.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param gene_id: (int) a sedDB gene id (same as a fasterDB gene id)
    :return: (string) the name of the exon
    """
    cursor = cnx.cursor()
    query = "SELECT DISTINCT gene_symbol FROM sed WHERE gene_id = ?"
    cursor.execute(query, (gene_id, ))
    res = cursor.fetchall()
    gene_name = [val[0] for val in res]
    return gene_name[0]


def get_splicing_factor_name(cnx):
    """
    Get the name of every splicing factor in splicing lore.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :return: (list of string) list of splicing factor name
    """
    cursor = cnx.cursor()
    query = "SELECT DISTINCT sf_name FROM rnaseq_projects"
    cursor.execute(query)
    res = cursor.fetchall()
    sf_name = [val[0] for val in res]
    return sf_name


def get_projects_links_to_a_splicing_factor(cnx, sf_name):
    """
    Get the id of every projects corresponding to a particular splicing factor.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_name: (string) splicing factor name
    :return: (list of int) list of id_project
    """
    cursor = cnx.cursor()
    query = "SELECT id FROM rnaseq_projects WHERE sf_name = ?"
    cursor.execute(query, (sf_name,))
    res = cursor.fetchall()
    idp = [val[0] for val in res if val[0] not in group_factor.bad_id_projects]
    return idp


def get_ase_events(cnx, id_project):
    """
    Get every exon regulated in a particular project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :return: (list of tuple of one str 2 int) each sublist corresponds to an exon (\
    exon_regulation + gene_id + exon_position on gene)
    """
    cursor = cnx.cursor()
    query = """SELECT delta_psi, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = ?
               AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
               AND pvalue_glm_cor <= 0.05"""
    cursor.execute(query, (id_project,))
    res = cursor.fetchall()
    if len(res) == 0:
            query = """SELECT delta_psi, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = ?
               AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
               AND pvalue <= 0.05"""
            cursor.execute(query, (id_project,))
            res = cursor.fetchall()
    nres = []
    for exon in res:
        nexon = list(exon[1:3])
        if exon[0] < 0:
            nexon = ["down"] + nexon
        else:
            nexon = ["up"] + nexon
        nres.append(nexon)
    return nres


def washing_events(exon_list):
    """
    Remove redundant exons or remove exons showing different regulation.

    :param exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + \
    exon_position on gene + exon_regulation). \
    Every exon regulated by a splicing factor in different projects.
    :return: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects without redundancy.
    """
    replace_dic = {"up": "down", "down": "up"}
    dic = {}
    prefix_list = []
    for exon in exon_list:
        exon_prefix = "%s_%s" % (exon[1], exon[2])
        exon_name = "%s_%s" % (exon[0], exon_prefix)
        if exon_name not in dic:
            if exon_prefix not in prefix_list:
                dic[exon_name] = 1
                prefix_list.append(exon_prefix)
            else:
                reverse_name = exon_name.replace(exon[0], replace_dic[exon[0]])
                if reverse_name in dic:
                    del(dic[reverse_name])
                # Else : the exon was deleted before because of a different regulation
        else:
            dic[exon_name] += 1
    # creation of the new list of exons
    new_exon_list = []
    for key in dic:
        my_exon = key.split("_")
        new_exon_list.append(my_exon)
    return new_exon_list


def get_every_events_4_a_sl(cnx, sf_name, regulation):
    """
    Get every splicing events for a give splicing factor.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_name: (string) the name of a splicing factor
    :param regulation: (string) up or down
    :return: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects.
    """
    exons_list = []
    id_projects = get_projects_links_to_a_splicing_factor(cnx, sf_name)
    for id_project in id_projects:
        ase_event = get_ase_events(cnx, id_project)
        exons_list += ase_event

    washed_exon_list = washing_events(exons_list)
    reg_exon_list = []
    for exon in washed_exon_list:
        if exon[0] == regulation:
            reg_exon_list.append(exon[1:])
    return reg_exon_list


def get_events_4_a_sl(cnx, sf_name, regulation):
    """
    Get every splicing events for a give splicing factor. Every exon down or up saw at least once will be reported \
    even if for a splicing factor some exons are up-and down regulated in differents project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_name: (string) the name of a splicing factor
    :param regulation: (string) up or down
    :return: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects.
    """
    exons_list = []
    id_projects = get_projects_links_to_a_splicing_factor(cnx, sf_name)
    for id_project in id_projects:
        ase_event = get_ase_events(cnx, id_project)
        exons_list += ase_event

    reg_exon_list = []
    for exon in exons_list:
        if exon[0] == regulation:
            reg_exon_list.append(exon[1:])
    # removing redundant exons
    washed_exon_list = washing_events_all(reg_exon_list)

    return washed_exon_list


def washing_events_all(exon_list):
    """
    Remove redundant exons

    :param exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon \
    (gene_id + exon_position on gene + exon_regulation). \
    Every exon regulated by a splicing factor in different projects.
    :return: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects without redundancy.
    """
    dic = {}
    for exon in exon_list:
        exon_name = "%s_%s" % (exon[0], exon[1])
        if exon_name not in dic:
            dic[exon_name] = 1
        else:
            dic[exon_name] += 1
    # creation of the new list of exons
    new_exon_list = []
    for key in dic:
        my_exon = key.split("_")
        new_exon_list.append(my_exon)
    return new_exon_list


def get_exon_regulated_by_sf(cnx, regulation):
    """
    Get the exons ``regulation`` regulated by a splicing factors.

    :param cnx: (sqlite3 connect object) connection to sed database
    :param regulation: (string) up or down
    :return: (list of list of 2 int) list of exons regulated by a splicing factor
    """
    name_projects = group_factor.get_wanted_sf_name(None)
    exon_list = []
    for sf_name in name_projects:
        exon_list += get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_list = np.unique(exon_list, axis=0).tolist()
    exon_list = [list(map(int, exon)) for exon in exon_list]
    print("Number of exons regulated by a splicing factor : %s" % len(exon_list))
    return exon_list
