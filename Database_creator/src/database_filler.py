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
import database_creator



# Functions
def connection():
    """
    :return: (pymysql object to connected to the fasterDB database)
    """
    cnx = pymysql.connect(user=conf.user, password=conf.password, host=conf.host, database=conf.fasterDB)
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
    cursor.executemany("INSERT INTO introns VALUES (?, ?, ?, ?)", result)
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
    SELECT id_gene, est_site_donor, exon_pos, `force`
    FROM force_splicing_site WHERE est_alternatif=0
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor = new_db.cursor()
    cursor.executemany("INSERT INTO force_splicing_site VALUES (?, ?, ?, ?)", result)
    new_db.commit()


def fill_exon_genomiques_table(new_db):
    """
    Fill the table **exon_genomiques** in ``new_db``

    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = new_db.cursor()
    query = """
    SELECT t1.id_gene, t1.pos_on_gene, t1.start_on_gene, t1.end_on_gene, t1.exon_type, t1.cds_start_on_gene, t1.cds_end_on_gene, t1.offset_before_exon, t1.offset_after_exon, t1.force_donor, t2.force as force_acceptor
    FROM(
        SELECT t1.id_gene, t1.pos_on_gene, t1.start_on_gene, t1.end_on_gene, t1.exon_type, t1.cds_start_on_gene, t1.cds_end_on_gene, t1.offset_before_exon, t1.offset_after_exon, t2.force as force_donor
        FROM exon_partial t1 LEFT JOIN (SELECT * FROM force_splicing_site WHERE is_donor = 1) t2
        ON t1.id_gene = t2.id_gene
        AND t1.pos_on_gene = t2.pos_on_gene) t1
    LEFT JOIN (SELECT * FROM force_splicing_site WHERE is_donor = 0) t2
    ON t1.id_gene = t2.id_gene
    AND t1.pos_on_gene = t2.pos_on_gene
    """
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.executemany("INSERT INTO exons VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", result)
    new_db.commit()


def remove_exon_patial_and_force_splicing_site(new_db):
    """
    This function will destroy the tables **force_splicing_site** and **exon_partial** tables
    :param new_db: (sqlite3 object) connection to ``new_db``
    """
    cursor = new_db.cursor()
    cursor.execute("DROP TABLE force_splicing_site;")
    cursor.execute("DROP TABLE exon_partial;")
    new_db.commit()


def database_maker():
    """
    :return:  Create the fasterDB lite database
    """
    print("database_creation")
    base_name = database_creator.database_creator()
    print("establishing connextion betwwen 2 cards")
    new_db = database_creator.new_db_connection(base_name)
    cnx = connection()
    print("filling genes content tables")
    fill_gene_table_content(cnx, new_db)
    print("filling intron table")
    fill_intron_table(cnx, new_db)
    print("filling partial exon table...")
    fill_exon_partial_table(cnx, new_db)
    print("force table")
    fill_force_table(cnx, new_db)
    print("filling full force table")
    fill_exon_genomiques_table(new_db)
    print("removing tables exon_partial and force_splicing_site")
    remove_exon_patial_and_force_splicing_site(new_db)
    print("succefully ending...")
    cnx.close()
    new_db.close()

if __name__ == "__main__":
    database_maker()
