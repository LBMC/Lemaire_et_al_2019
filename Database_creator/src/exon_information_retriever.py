#!/usr/bin/python3.5

# coding: utf-8

# set environment

import exon_class
import sqlite3
from database_creator import base_name, out_path
import sys
import re
import conf
import pymysql

# functions
def connection_sl():
    """
    :return: (pymysql object to connected to the splicing lore database)
    """
    cnx = pymysql.connect(user=conf.user, password=conf.password, host=conf.sl_host, database=conf.fasterDB)
    return cnx


def fasterdbl_connection(fasterdb):
    """
    :param fasterdb: (string) path to fasterDB database

    Allow connection to fasterDB Lite
    """
    return sqlite3.connect(fasterdb)


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
        if exon.upstream_exon.donor is not None and exon.upstream_exon.donor != 0 and exon.donor is not None:
            relative_donor_upstream = round(((exon.donor - exon.upstream_exon.donor) / exon.upstream_exon.donor) * 100, 1)
        else:
            relative_donor_upstream = None
        if exon.downstream_exon.donor is not None and exon.downstream_exon.donor != 0 and exon.donor is not None:
            relative_donor_downstream = round(((exon.donor - exon.downstream_exon.donor) / exon.downstream_exon.donor) * 100, 1)
        else:
            relative_donor_downstream = None

        if exon.upstream_exon.acceptor is not None and exon.upstream_exon.acceptor != 0 and exon.acceptor is not None:
            relative_acceptor_upstream = round(((exon.acceptor - exon.upstream_exon.acceptor) / exon.upstream_exon.acceptor) * 100, 1)
        else:
            relative_acceptor_upstream = None
        if exon.downstream_exon.acceptor is not None and exon.downstream_exon.acceptor != 0 and exon.acceptor is not None:
            relative_acceptor_downstream = round(((exon.acceptor - exon.downstream_exon.acceptor) / exon.downstream_exon.acceptor) * 100, 1)
        else:
            relative_acceptor_downstream = None
        cur_list = [exon.gene.name, exon.gene.id, exon.position, exon.type, exon.gene.length,
                    exon.gene.nb_intron, exon.gene.median_intron_size, exon.gene.iupac, exon.gene.dnt,
                    exon.upstream_exon.length, exon.length, exon.downstream_exon.length,
                    exon.upstream_intron.length, exon.downstream_intron.length,
                    exon.upstream_exon.acceptor, exon.acceptor, exon.downstream_exon.acceptor,
                    exon.upstream_exon.donor, exon.donor, exon.downstream_exon.donor, exon.iupac, exon.dnt,
                    exon.upstream_intron.iupac_dist, exon.upstream_intron.dnt_dist,
                    exon.upstream_intron.iupac_proxi, exon.upstream_intron.dnt_proxi,
                    exon.downstream_intron.iupac_proxi, exon.downstream_intron.dnt_proxi,
                    exon.downstream_intron.iupac_dist, exon.downstream_intron.dnt_dist,
                    relative_donor_upstream, relative_donor_downstream, relative_acceptor_upstream,
                    relative_acceptor_downstream]
        list_tuple.append(cur_list)
    return list_tuple


def sed_connection(seddb):
    """
    :param seddb: (string) path to sed database

    Connection to SED database
    """
    return sqlite3.connect(seddb)


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
               iupac_gene VARCHAR(50),
               dnt_gene VARCHAR(80),
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
               dnt_exon VARCHAR(80),
               iupac_upstream_intron_dist VARCHAR(50),
               dnt_upstream_intron_dist VARCHAR(80),
               iupac_upstream_intron_proxi VARCHAR(50),
               dnt_upstream_intron_proxi VARCHAR(80),
               iupac_downstream_intron_proxi VARCHAR(50),
               dnt_downstream_intron_proxi VARCHAR(80),
               iupac_downstream_intron_dist VARCHAR(50),
               dnt_downstream_intron_dist VARCHAR(80),
               relative_donor_upstream FLOAT,
               relative_donor_downstream FLOAT,
               relative_acceptor_upstream FLOAT,
               relative_acceptor_downstream FLOAT,
               PRIMARY KEY(gene_id, exon_pos));
               """
    cursor.execute(query)
    sed_cnx.commit()


def creation_rnaseq_projects_table(sed_cnx):
    """
    Create a rnaseq_projects table in ``sed_cnx``

    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    """
    cursor = sed_cnx.cursor()
    query = """
    CREATE TABLE rnaseq_projects (
        id INT,
        project_name VARCHAR(50) NOT NULL,
        source_db VARCHAR(10) NOT NULL,
        db_id_project VARCHAR(20) NOT NULL,
        sf_name VARCHAR(25) NOT NULL,
        cl_name VARCHAR(25) NOT NULL,
        PRIMARY KEY (id)
    );
    """
    cursor.execute(query)
    sed_cnx.commit()


def creation_ase_event_tmp_table(sed_cnx):
    """
    Create an ase_event_tmp table in ``sed_cnx``

    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    """
    cursor = sed_cnx.cursor()
    query = """
    CREATE TABLE ase_event_tmp (
        id INT,
        id_project INT NOT NULL,
        gene_symbol VARCHAR(20) NOT NULL,
        exon_skipped INT NOT NULL,
        chromosome VARCHAR(1) NOT NULL,
        start INT NOT NULL,
        stop INT NOT NULL,
        exons_flanquants VARCHAR(15),
        delta_psi FLOAT,
        pvalue FLOAT,
        pvalue_glm_cor FLOAT,
        PRIMARY KEY (id)
    );
    """
    cursor.execute(query)
    sed_cnx.commit()


def creation_ase_event_table(sed_cnx):
    """
    Create an ase_event table in ``sed_cnx``

    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    """
    cursor = sed_cnx.cursor()
    query = """
    CREATE TABLE ase_event (
        id INT,
        id_project INT NOT NULL,
        gene_id INT NOT NULL,
        gene_symbol VARCHAR(20) NOT NULL,
        exon_skipped INT NOT NULL,
        chromosome VARCHAR(1) NOT NULL,
        start INT NOT NULL,
        stop INT NOT NULL,
        exons_flanquants VARCHAR(15),
        delta_psi FLOAT,
        pvalue FLOAT,
        pvalue_glm_cor FLOAT,
        PRIMARY KEY (id),
        FOREIGN KEY (id_project) REFERENCES rnaseq_projects(id),
        FOREIGN KEY (gene_id, exon_skipped) REFERENCES sed(gene_id, exon_pos)
    );
    """
    cursor.execute(query)
    sed_cnx.commit()


def sed_filler(sed_cnx, list_tuple):
    """
    Fill the **sed table**

    :param sed_cnx: (sqlite3 object) contains every information we need to connect to \
    SED database.
    :param list_tuple: (list of **list info exon**). **list info exon** contains every information \
    necessary for an exon.
    """
    cursor = sed_cnx.cursor()
    cursor.executemany(
        "INSERT INTO sed VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        list_tuple)
    sed_cnx.commit()
    # creation of an index on gene_symbol and exon_pos
    query = """ CREATE INDEX sed_index on sed(gene_symbol, exon_pos);"""
    cursor.execute(query)
    sed_cnx.commit()


def fill_rnaseq_projects_content(cnx, sed_cnx):
    """
    Fill the table **rnaseq_projects** in ``sed database``.

    :param cnx: (pymysql object) connection to splicing lore
    :param sed_cnx: (sqlite3 object) connection to ``sed database``
    """
    cursor = cnx.cursor()
    query = """
    SELECT t1.id, t1.project_name, t1.source_db, t1.db_id_project, t2.name, t3.name_CL
    FROM rnaseq_projects_SF t1, Splicing_factors t2, Cell_line t3
    WHERE t1.id_SF = t2.id_SF
    AND t1.id_cell_line = t3.id_cell_line
    AND t1.show_in_website = 1
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = sed_cnx.cursor()
    cursor.executemany("INSERT INTO rnaseq_projects VALUES (?, ?, ?, ?, ?, ?)", result)
    sed_cnx.commit()


def fill_ase_event_tmp_content(cnx, sed_cnx):
    """
    Fill the table **ase_event_tmp** in ``sed database``.

    :param cnx: (pymysql object) connection to splicing lore
    :param sed_cnx: (sqlite3 object) connection to ``sed database``
    """
    cursor = cnx.cursor()
    # selection of ASE events show in splicing Lore.
    query = """
    SELECT t1.id, t1.id_project, t1.gene_symbol, t1.exon_skipped, t1.coordonnees,
    t1.exons_flaquants, t1.DeltaPSI, t1.pvalue, t1.pvalue_glm_corr
    FROM Events_ASE t1, rnaseq_projects_SF t2
    WHERE t1.id_project = t2.id
    AND t2.show_in_website = 1
    """
    cursor.execute(query)
    result = cursor.fetchall()
    new_result = []
    for row in result:
        new_row = [row[i] for i in range(4)]
        coord = re.split(":|-", row[4])
        new_row.append(coord[0])
        new_row.append(int(coord[1]))
        new_row.append(int(coord[2]))
        new_row += [row[i] for i in range(5, len(row), 1)]
        new_result.append(new_row)

    cursor = sed_cnx.cursor()
    my_tuple = str(tuple((["?"] * len(new_result[0])))).replace("'", "")
    cursor.executemany("INSERT INTO ase_event_tmp VALUES %s" % my_tuple, new_result)
    sed_cnx.commit()


def fill_ase_event_content(sed_cnx, seddb, fasterdblite):
    """
    Fill the table **ase_event** in ``sed_cnx``.

    :param sed_cnx: (sqlite3 object) connection to ``sed database``
    :param seddb: (string) path to sed database
    :param fasterdblite: (string) path to fasterdblite database
    """
    query = """
    SELECT t1.id, t1.id_project, t2.id_gene, t1.gene_symbol,
    t1.exon_skipped, t1.chromosome, t1.start, t1.stop, t1.exons_flanquants, t1.delta_psi, t1.pvalue, t1.pvalue_glm_cor
    FROM sed.ase_event_tmp t1, fasterdb.exons t2, fasterdb.genes t3
    WHERE t2.id_gene = t3.id
    AND t1.chromosome = t2.chromosome
    AND t1.start = t2.start_on_chromosome
    AND t1.stop = t2.end_on_chromosome
    AND t3.official_symbol = t1.gene_symbol
    AND t1.exon_skipped = t2.pos_on_gene
    """

    cursor = sed_cnx.cursor()
    cursor.execute("ATTACH DATABASE ? as sed", (seddb,))
    cursor.execute("ATTACH DATABASE ? as fasterdb", (fasterdblite,))
    cursor.execute(query)
    result = cursor.fetchall()
    my_tuple = str(tuple((["?"] * len(result[0])))).replace("'", "")
    cursor.executemany("INSERT INTO ase_event VALUES %s" % my_tuple, result)
    sed_cnx.commit()
    query = """ CREATE INDEX sed_project on ase_event(id_project);"""
    cursor.execute(query)
    sed_cnx.commit()
    query = """ CREATE INDEX sed_regulation on ase_event(id_project, delta_psi, pvalue);"""
    cursor.execute(query)
    sed_cnx.commit()
    query = """ CREATE INDEX sed_regulation_glm on ase_event(id_project, delta_psi, pvalue_glm_cor);"""
    cursor.execute(query)
    sed_cnx.commit()
    query = """ CREATE INDEX project on rnaseq_projects(sf_name, cl_name);"""
    cursor.execute(query)
    sed_cnx.commit()


def remove_ase_event_tmp(sed_cnx):
    """
    Destroy the tables **ase_event_tmp** table.

    :param sed_cnx: (sqlite3 object) connection to ``sed database``
    """
    cursor = sed_cnx.cursor()
    cursor.execute("DROP TABLE ase_event_tmp;")
    sed_cnx.commit()


def main():
    # debug mode
    debug = 0  # 1 = enabled , 0 disabled
    fasterdblite = out_path + base_name
    seddb = out_path + "sed.db"
    cnx = fasterdbl_connection(out_path + base_name)
    info_list = exon_finder(cnx)
    exon_list = get_exon_info(cnx, info_list, debug)
    list_tuple = get_exon_tuple(exon_list)
    cnx.close()
    cnx = connection_sl()
    sed_cnx = sed_connection(seddb)
    print("Creation of SED table")
    create_sed_exon_table(sed_cnx)
    print("Creation of rnaseq_projects table")
    creation_rnaseq_projects_table(sed_cnx)
    print("Creation of ase_event_tmp table")
    creation_ase_event_tmp_table(sed_cnx)
    print("Creation of ase_event table")
    creation_ase_event_table(sed_cnx)
    print("Filling sed table")
    sed_filler(sed_cnx, list_tuple)
    print("Filling rnaseq_projects table")
    fill_rnaseq_projects_content(cnx, sed_cnx)
    print("Filling ase_event_tmp table")
    fill_ase_event_tmp_content(cnx, sed_cnx)
    print("Filling ase_event table")
    fill_ase_event_content(sed_cnx, seddb, fasterdblite)
    print("Removing ase_event_tmp table")
    remove_ase_event_tmp(sed_cnx)
    print("closing connections")
    cnx.close()
    sed_cnx.close()
    print("successfully ended ! ")


if __name__ == "__main__":
    main()
