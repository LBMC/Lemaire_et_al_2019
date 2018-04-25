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
    Fill the table gene in ``new_db``

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





