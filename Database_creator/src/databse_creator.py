#!/bin/usr/python3.5

"""
Description:

    This script will be use to create a database containing 3 tables : \

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
    cnx = pymysql.connect(user=conf.user, password=conf.password, host=conf.host, database="Nicolas")
    return cnx


def new_db_connection(path2db):
    """
    :param path2db: (string) path where the database will be or is located
    """
    return sqlite3.connect(path2db + "fasterDB_lite.db")



def creation_of_gene_table(new_db):
    """
    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE genes (
        id smallint(5) NOT NULL,
        official_symbol VARCHAR(17) NOT NULL,
        chromosome VARCHAR(2) NOT NULL,
        strand tinyint(1) NOT NULL,
        start_on_chromosome int(10) NOT NULL,
        end_on_chromosome int(10) NOT NULL,
        sequence mediumtext NOT NULL,
        PRYMARY KEY (id)
    );
    """
    cursor.execute(query)
    new_db.commit()
