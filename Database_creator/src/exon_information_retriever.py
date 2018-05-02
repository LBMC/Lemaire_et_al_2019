#!/usr/bin/python3.5

# coding: utf-8

# set environment

import exon_class
import sqlite3
from database_creator import base_name, out_path
import sys

# functions
def fasterdbl_connection():
    """
    Allow connection to fasterDB Lite
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
    return cursor.fetchall()


def get_exon_info(cnx, info_list, debug):
    """
    Get every information we need on an exon

    :param debug: (int) 0 no debug, 1 debug mode
    :param cnx: (sqlite3 object) return all the information we need to connect to FasterDB lite
    :param info_list: (list of list of string and int and int) each sublist contains \
    a string : gene_symbol and 2 int : the gene_id and the exobn position on gene respectively
    :return: (a list of ExonClass object) list of exons
    """
    print("Getting exons information !")
    exon_list = []
    exon_class.set_debug(debug)
    count = 0
    ll = str(len(info_list))
    for exon_info in info_list:
        exon_list.append(exon_class.ExonClass(cnx, exon_info[0], exon_info[1], exon_info[2]))
        count += 1
        percent = round(float(count) / len(info_list) * 100, 1)
        sys.stdout.write("Progression : " + str(count) + " / " + ll + " - " + str(percent) + " %\r")
        sys.stdout.flush()
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


def sed_connection():
    """
    Connection to SED database
    """
    return sqlite3.connect(out_path + "sed.db")


def create_sed_exon_table(sed_cnx):
    """
    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    """
    cursor = sed_cnx.cursor()
    query = """CREATE TABLE sed (
               gene_symbol VARCHAR(17) NOT NULL,
               gene_id INT(10) NOT NULL,
               exon_pos INT(10) NOT NULL,
               exon_type VARCHAR(10),
               gene_size INT(10) NOT NULL,
               nb_intron_gene INT NOT NULL,
               median_intron_size INT,
               iupac_gene VARCHAR(50) NOT NULL,
               upstream_exon_size INT,
               exon_size INT,
               downstream_exon_size INT,
               upstream_intron_size INT,
               downstream_intron_size INT,
               force_acceptor_upstream_exon FLOAT,
               force_acceptor FLOAT,
               force_acceptor_downstream_exon FLOAT,
               force_donor_upstream_exon FLOAT,
               force_donor FLOAT,
               force_donor_downstream_exon FLOAT,
               iupac_exon VARCHAR(50),
               iupac_upstream_intron_dist VARCHAR(50),
               iupac_upstream_intron_proxi VARCHAR(50),
               iupac_downstream_intron_proxi VARCHAR(50),
               iupac_downstream_intron_dist VARCHAR(50),
               PRIMARY KEY(gene_id, exon_pos));
               """
    cursor.execute(query)
    sed_cnx.commit()


def sed_filler(sed_cnx, list_tuple):
    """
    Fill the **sed database**

    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    :param list_tuple: (list of **list info exon**). **list info exon** contains every information \
    necessary for an exon.
    """
    cursor = sed_cnx.cursor()
    cursor.executemany(
        "INSERT INTO sed VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        list_tuple)
    sed_cnx.commit()


def main():
    # debug mode
    debug = 0  # 1 = enabled , 0 disabled

    cnx = fasterdbl_connection()
    info_list = exon_finder(cnx)
    exon_list = get_exon_info(cnx, info_list, debug)
    list_tuple = get_exon_tuple(exon_list)
    sed_cnx = sed_connection()
    create_sed_exon_table(sed_cnx)
    sed_filler(sed_cnx, list_tuple)
    cnx.close()
    sed_cnx.close()


if __name__ == "__main__":
    main()
