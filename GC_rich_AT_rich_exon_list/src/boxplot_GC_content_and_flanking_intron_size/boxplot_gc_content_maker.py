#!/usr/bin/env python3

# -*- coding utf-8 -*-

"""
Description: This program will a graphics with the gc content of AT and GC group of exons (same \
for the gene containing those exons)
"""

import numpy as np


def get_exon_control_gc_content(cnx, exon_type, exon2remove):
    """
    Get the GC content of every exon with ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sedb database
    :param exon_type: (string) the ``exon_type`` of interest
    :param exon2remove: (list of list of 2 int) list of exons 2 remove
    :return: (list of float) list of the gc content of every exons with the ``exon_type`` of interest
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos, iupac_exon
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos, iupac_exon
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    print("Control exons before removing bad ones : %s" % len(tuple_list))
    gc_content = [exon[2].split(";")[4] for exon in tuple_list if [exon[0], exon[1]] not in exon2remove]
    print("Control exons after removing those regulated by a sf : %s" % len(gc_content))
    return gc_content


def get_gene_control_gc_content(cnx, exon_type, gene2remove):
    """
    Get the GC content of every gene containing exons with ``exon_type``
    :param cnx: (sqlite3 connect object) connection to sedb database
    :param exon_type: (string) the ``exon_type`` of interest
    :param gene2remove: (string) gene regulated by a splicing factor
    :return: (list of float) list of the gc content of every exons with the ``exon_type`` of interest
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT DISTINCT iupac_gene, gene_id
                       FROM sed
                       WHERE exon_type LIKE '%{}%'
                       """.format(exon_type)
    else:
        query = """SELECT DISTINCT iupac_gene, gene_id
                       FROM exons
                    """
    cursor.execute(query)
    tuple_list = cursor.fetchall()
    gc_content = []
    for iupac in tuple_list:
        if iupac[1] not in gene2remove:
            if iupac[0]:
                gc_content.append(iupac[0].split(";")[4])
            else:
                gc_content.append(None)
    print("NB gene not containing sf-regulated exons : %s" % len(gc_content))
    # turn tuple into list
    return gc_content


def calculate_exon_gc_content(cnx, gene_id, exon_pos):
    """
    Get the GC content of an exon having the following ``gene_id`` and ``exon_pos``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param gene_id: (int) the fasterdb ``gene_id`` of a gene
    :param exon_pos: (int) the fasterDb position of an exon within the gene``gene_id``
    :return: (float) the gc content of the exon identified by ``gene_id`` and ``exon_position``
    """
    cursor = cnx.cursor()
    query = "SELECT iupac_exon FROM sed WHERE gene_id = ? and exon_pos = ?"
    cursor.execute(query, (gene_id, exon_pos,))
    res = cursor.fetchone()
    return res[0].split(";")[4]


def calculate_gene_gc_content(cnx, gene_id):
    """
    Get the GC content of a gene containing with the following ``gene_id``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param gene_id: (int) the fasterdb ``gene_id`` of a gene
    :return: (float) the gc content of the gene with the ``gene_id``
    """
    cursor = cnx.cursor()
    query = "SELECT DISTINCT iupac_gene FROM sed WHERE gene_id = ?"
    cursor.execute(query, (gene_id,))
    res = cursor.fetchall()
    if len(res) > 1:
        print("WARNING : the gene selected has different gc content possible : %s" % str(res))
    if res[0][0]:
        return res[0][0].split(";")[4]
    else:
        return None


def extract_exon_gc_content_from_file(cnx, filename):
    """
    Extract the gc_ content of a list of exon within a file ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) the name of the file containing  exons
    :return: (list of float) the list of the gc content of the exon in ``filename``
    """
    list_gc = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            list_gc.append(calculate_exon_gc_content(cnx, line[0], line[1]))
            line = in_file.readline()
    return list_gc


def extract_gene_gc_content_from_file(cnx, filename, gene2remove):
    """
    Extract the gc content of genes of a list of exon within a file ``filename``
    :param cnx: (sqlite3 connect object) connection to sed database
    :param filename: (string) the name of the file containing  exons
    :param gene2remove: (list of string) list of gene 2 remove
    :return: (list of float) the list of the gc content of the gene in ``filename``
    """
    list_gene = []
    with open(filename, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.split("\t")
            list_gene.append(line[0])
            line = in_file.readline()
    list_gene = list(np.unique(list_gene))
    print("list gene : %s" % len(list_gene))
    list_gene = [g for g in list_gene if g not in gene2remove]
    print("number of gc gene after removing common with at and GC exons : %s" % len(list_gene))
    list_gc = extract_gene_gc_content_from_list(cnx, list_gene)
    return list_gc


def extract_exon_gc_content_from_list(cnx, exon_list):
    """
    Get the gc content from a list of exons
    :param cnx: (sqlite3 connection object) connection to sed database
    :param exon_list: (list of 2 int) list of exons identified by its gene_id and position wihtin this gene
    :return: (list of float) the list of gc content of every exons within
    """
    list_gc = []
    for exon in exon_list:
        list_gc.append(calculate_exon_gc_content(cnx, exon[0], exon[1]))
    return list_gc


def extract_gene_gc_content_from_list(cnx, gene_list):
    """
    Get the gc content from a list of gene
    :param cnx: (sqlite3 connection object) connection to sed database
    :param gene_list: (list of 1 int) list of gene_id
    :return: (list of float) the list of gc content of every gene within ``gene_list``
    """
    list_gc = []
    for gene in gene_list:
        val = calculate_gene_gc_content(cnx, gene)
        if val:
            list_gc.append(val)
        else:
            list_gc.append(None)
    return list_gc