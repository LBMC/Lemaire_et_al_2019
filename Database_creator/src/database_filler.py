#!/bin/usr/python3.5

"""
Description:

    This script will be use to create and fille a database containing 3 empty tables : \

    * exon_genomiques
    * intron_genomiques
    * genes
    Those data are taken from fasterDB
"""

# sets the environment
import pymysql
import sqlite3
import conf



# Functions
def connection(dest):
    """
    :param dest: (string) the database for which we want a connection
    :return: (pymysql object to connected to the wanted database)
    """
    if dest == "fasterDB":
        db = conf.fasterDB
    else:
        db = conf.my_db
    cnx = pymysql.connect(user=conf.user, password=conf.password, host=conf.host, database=db)
    return cnx


def fill_gene_table_content(cnx, new_db):
    """
    Fill the table **genes** in ``new_db``

    :param cnx: (pymysql object) connection to fasterDB human
    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = cnx.cursor()
    query = """
    SELECT id, official_symbol, chromosome,	strand, start_sur_chromosome, end_sur_chromosome, sequence
    FROM genes;
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = new_db.cursor()
    cursor.executemany("INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?)", result)
    new_db.commit()


def fill_intron_table(cnx, new_db):
    """
    Fill the table **intron** in ``new_db``

    :param cnx: (pymysql object) connection to fasterDB human
    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = cnx.cursor()
    query = """
    SELECT id_gene, pos_sur_gene, start_sur_gene, end_sur_gene
    FROM introns_genomiques;
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = new_db.cursor()
    cursor.executemany("INSERT INTO intron_genomiques VALUES (?, ?, ?, ?)", result)
    new_db.commit()


def fill_exon_partial_table(cnx, new_db):
    """
    Fill the table **exon_partial** in ``new_db``

    :param cnx: (pymysql object) connection to fasterDB human
    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = cnx.cursor()
    query = """
    SELECT t1.id_gene, t1.pos_sur_gene, t1.start_sur_gene, t1.end_sur_gene, t1.exon_types, t2.fragment_start_on_gene, t2.fragment_end_on_gene, t2.offset_before_exon, t2.offset_after_exon
    FROM (
        SELECT t1.id_gene, t1.pos_sur_gene, t1.start_sur_gene, t1.end_sur_gene, t2.exon_types
        FROM fasterdb_humain.exons_genomiques t1, fasterdb_protein.hsapiens_exonsstatus_improved t2
        WHERE t1.id_gene = t2.id_gene
        AND t1.pos_sur_gene = t2.pos_sur_gene
    )t1 LEFT JOIN Nicolas.hsapiens_exonpeptides_filtered t2 ON t1.id_gene = t2.gene_id
    AND t1.pos_sur_gene = t2.exon_position_on_gene;
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = new_db.cursor()
    cursor.executemany("INSERT INTO exon_partial VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", result)
    new_db.commit()


def fill_force_table(cnx, new_db):
    """
    Fill the table **force_splicing_site** in ``new_db``

    :param cnx: (pymysql object) connection to fasterDB human
    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = cnx.cursor()
    query = """
    SELECT id_gene, est_site_donor, exon_pos, force, est_alternatif
    FROM force_splicing_site;
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = new_db.cursor()
    cursor.executemany("INSERT INTO force_splicing_site VALUES (?, ?, ?, ?, ?)", result)
    new_db.commit()
