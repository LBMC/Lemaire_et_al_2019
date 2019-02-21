#!/usr/bin/python3

# -*- coding utf-8 -*-


"""
Description: This program will a graphics with the flanking intron size of AT, GC CCE group of exons and the \
median_intron_size of the AT, GC, CCE genes (gene containing at least one exons of the AT, GC or CCE group)
"""

import numpy as np


def get_exon_control_min_flanking_intron_size(cnx, exon_type, exon2remove):
    """
    Get the min flanking intron size of every exons with the exon type ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_type: (string) the type of exon for which we want to get
    :param exon2remove: (string) the list of exon we want to remove
    :return: (list of float) the flanking intron size of exons with the exon type ``exon_tyep``
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos, upstream_intron_size, downstream_intron_size
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos, upstream_intron_size, downstream_intron_size
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    median_intron_size = [np.nanmin(np.array([size[2], size[3]], dtype=float)) for size in tuple_list
                          if [size[0], size[1]] not in exon2remove]
    print("Control exons after removing those regulated (down or up) by a splicing factor : %s" % len(median_intron_size))
    # turn tuple into list
    return median_intron_size


def get_gene_control_median_flanking_intron_size(cnx, exon_type, gene2remove):
    """
    Get the median intron size of every gene containing at least one exon type ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_type: (string) the type of exon for which we want to
    :param gene2remove: (list of int) list of genes to remove
    :return: (list of float) the medain intron size of exons with the exon type ``exon_type``
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT DISTINCT median_intron_size, gene_id
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT DISTINCT median_intron_size, gene_id
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    median_intron_size = [size[0] for size in tuple_list if size[0] and size[1] not in gene2remove]
    print("NB gene control median_intron_size retrieved (not regulated by SF) : %s" % len(median_intron_size))
    return median_intron_size


def calculate_exon_min_flanking_intron_size(cnx, gene_id, exon_pos):
    """
    Get the min flanking intron size of the exons with ``gene_id`` and ``exon_pos``
    :param cnx: (sqlite3 connect object) connection to sed database.
    :param gene_id: (int) the id of the gene containing the exons.
    :param exon_pos: (int) the id of the exon of interest
    :return: (float) the median flanking intron size of the exon with ``gene_id`` and ``exon_pos``
    """
    cursor = cnx.cursor()
    query = "SELECT upstream_intron_size, downstream_intron_size FROM sed WHERE gene_id = ? and exon_pos = ?"
    cursor.execute(query, (gene_id, exon_pos,))
    res = cursor.fetchone()
    return np.nanmin(np.array([res[0], res[1]], dtype=float))


def calculate_gene_median_intron_size(cnx, gene_id):
    """
    Get the median flanking intron size of the gene ``gene_id``
    :param cnx: (sqlite3 connect object) connection to sed database.
    :param gene_id: (int) the id of the gene
    :return: (float) the median flanking intron size of the gene ``gene_id``
    """
    cursor = cnx.cursor()
    query = "SELECT median_intron_size FROM sed WHERE gene_id = ?"
    cursor.execute(query, (gene_id,))
    res = cursor.fetchone()
    return res[0]


def extract_exon_min_flanking_intron_size_from_list(cnx, exon_list):
    """
    Get the min flanking intron size of every exon located in exons list
    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_list: (list of 2 int) list of exons identified by their gene_id and their position in the gene.
    :return: (list of float) the list of median flaking intron size
    """
    list_gc = []
    for exon in exon_list:
        list_gc.append(calculate_exon_min_flanking_intron_size(cnx, exon[0], exon[1]))
    return list_gc


def extract_gene_median_intron_size_from_list(cnx, gene_list):
    """
    Get the median flanking intron size of every exon located in exons list
    :param cnx: (sqlite3 connect object) connection to sed database
    :param gene_list: (list of int) list of gene_id
    :return: (list of float) the list of median intron size of every gene in ``exon_list``
    """
    list_gc = []
    for gene in gene_list:
        val = calculate_gene_median_intron_size(cnx, gene)
        if val:
            list_gc.append(val)
    return list_gc


def extract_exon_min_flanking_intron_size_from_file(cnx, filename):
    """
    Get the min flanking intron size of the exon within the file ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) the file where the exons are stored
    :return: (list of float) list of median intron size of the exon in ``filename``
    """
    median_intron_size = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            median_intron_size.append(calculate_exon_min_flanking_intron_size(cnx, line[0], line[1]))
            line = in_file.readline()
    return median_intron_size


def extract_gene_median_intron_size_from_file(cnx, filename):
    """
    Get the median intron size of the gene in ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) the file where the gene are stored
    :return: (list of float) list of median intron size of the exon in ``filename``
    """
    list_gene = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            list_gene.append(line[0])
            line = in_file.readline()
    list_gene = list(np.unique(list_gene))
    list_gc = extract_gene_median_intron_size_from_list(cnx, list_gene)
    return list_gc

