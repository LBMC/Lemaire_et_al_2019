#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is to find the gene size of every gene in a givene exon list
"""


import numpy as np



def get_control_gene_size(cnx, exon_type, gene2remove):
    """
    Get thegene_size of every gene containing at least one exon type ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param exon_type: (string) the type of exon for which we want to
    :param gene2remove: (string) the gene id we want to remove.
    :return: (list of float) the gene size of exons with the exon type ``exon_type``
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT DISTINCT gene_size, gene_id
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT DISTINCT gene_size, gene_id
                   FROM sed
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    gene_size = [size[0] for size in tuple_list if size[0] and size[1] not in gene2remove]
    print("NB control gene_size retrieved (not regulated by splicing factors): %s" % len(gene_size))
    # turn tuple into list
    return gene_size


def calculate_gene_size(cnx, gene_id):
    """
    Get the gene size of the gene ``gene_id``
    :param cnx: (sqlite3 connect object) connection to sed database.
    :param gene_id: (int) the id of the gene
    :return: (float) the median flanking intron size of the gene ``gene_id``
    """
    cursor = cnx.cursor()
    query = "SELECT DISTINCT gene_size FROM sed WHERE gene_id = ?"
    cursor.execute(query, (gene_id,))
    res = cursor.fetchone()
    return res[0]


def extract_gene_size_from_list(cnx, gene_list):
    """
    Get the gene size of every exon located in exons list
    :param cnx: (sqlite3 connect object) connection to sed database
    :param gene_list: (list of int) list of gene_id
    :return: (list of float) the list of gene size of every gene in ``exon_list``
    """
    list_gc = []
    for gene in gene_list:
        val = calculate_gene_size(cnx, gene)
        if val:
            list_gc.append(val)
    return list_gc


def extract_gene_size_from_file(cnx, filename):
    """
    Get the gene size of the gene in ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) the file where the gene are stored
    :return: (list of float) list of gene size of the exon in ``filename``
    """
    list_gene = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            list_gene.append(line[0])
            line = in_file.readline()
    list_gene = list(np.unique(list_gene))
    list_gc = extract_gene_size_from_list(cnx, list_gene)
    return list_gc

