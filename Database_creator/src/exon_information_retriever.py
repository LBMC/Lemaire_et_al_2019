#!/usr/bin/python3.5

# coding: utf-8

# set environment

import exon_class
import sqlite3
import os
from database_creator import base_name, out_path
# debug mode
debug = 1 # 1 = enabled , 0 disabled


# functions
def fasterdbl_connection():
    """
    :param file_name: (string) file where the database will be created
    """
    return sqlite3.connect(out_path + base_name)


def exon_finder(cnx):
    """
    Find for every exons in fasterDB lite, its official gene symbol, its \
    gene id and its position within  a gene.

    :param cnx: (sqlite3 object) allows connection to fasterDB lite
    :return: (list of tuple of string int int). Each sublist of the list \
    contains:
        * the gene symbol of the gene that contains the exon
        * the gene id of the gene that contains the exon
        * the position of tne exon within the gene
    """
    cursor = cnx.cursor()
    query = """SELECT t2.official_symbol , t1.id_gene, t1.pos_on_gene
               FROM genes t2, exons t1
               WHERE t2.id = t1.id_gene"""
    cursor.execute(query)
    return(cursor.fetchall()[0:3])


def get_exon_info(cnx, info_list, debug):
    """
    Get every information we need on an exon

    :param debug: (int) 0 no debug, 1 debug mode
    :param cnx: (sqlite3 object) return all the information we need to connect to FasterDB lite
    :param info_list: (list of list of string and int and int) each sublist contains \
    a string : gene_symbol and 2 int : the gene_id and the exobn position on gene respectively
    :return: (a list of ExonClass object) list of exons
    """
    exon_list = []
    exon_class.set_debug(debug)
    for exon_info in info_list:
        exon_list.append(exon_class.ExonClass(cnx, exon_info[0], exon_info[1], exon_info[2]))
    return exon_list


def get_exon_tuple(exon_list):
    """

    :param exon_list: (a list of ExonClass object) list of exons
    :return: list of **list info exon**. **list info exon** contains every information \
    necessary for an exon.
    """
    list_tuple = []
    for exon in exon_list:
        cur_list = [exon.gene.name, exon.gene.id, exon.position, exon.type, exon.gene.length,
                    exon.gene.nb_intron, exon.gene.median_intron_size, exon.gene.iupac,
                    exon.upstream_exon.length, exon.length, exon.downstream_exon.length,
                    exon.upstream_intron.length, exon.downstream_intron.length,
                    exon.upstream_exon.acceptor, exon.acceptor, exon.downstream_exon.acceptor,
                    exon.upstream_exon.donor, exon.donor, exon.downstream_exon.donor, exon.iupac,
                    exon.upstream_intron.iupac_dist, exon.upstream_intron.iupac_proxi,
                    exon.downstream_intron.iupac_proxi, exon.downstream_intron.iupac_dist]
        list_tuple.append(cur_list)
    return list_tuple


def write_result(list_tuple):
    """

    :param list_tuple: (list of **list info exon**). **list info exon** contains every information \
    necessary for an exon.
    """
    with open(out_path + "info_exon.txt", "w") as outfile:
        header = "gene_symbol, gene_id, exon_pos, exon_type, gene_size, nb_intron_gene, median_intron_size, " \
                 "iupac_gene, upstream_exon_size, exon_size, downstream_eoxn_size, upstream_intron_size, " \
                 "downstream_intron_size, upstream_exon_acceptor, exon_acceptor, downstream_exon_acceptor, " \
                 "upstream_exon_donor, exon_donor, downstream_exon_donor, iupac_exon, iupac_upstream_intron_dist, " \
                 "iupac_upstream_intron_proxi, iupac_downstream_intron_proxi, iupac_downstream_intron_dist\n"
        outfile.write(header.replace(", ", "\t"))
        for mytuple in list_tuple:
            outfile.write(str(mytuple).replace("[", "").replace("]", "").replace(", ", "\t") + "\n")



def main():
    cnx = fasterdbl_connection()
    info_list = exon_finder(cnx)
    exon_list = get_exon_info(cnx, info_list, debug)
    list_tuple = get_exon_tuple(exon_list)
    write_result(list_tuple)






