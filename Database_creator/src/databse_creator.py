#!/bin/usr/python3.5

"""
Description:

    This script will be use to create a database containing 3 empty tables : \

    * exon_genomiques
    * intron_genomiques
    * genes
    Those data are taken from fasterDB
"""

# sets the environment
import pymysql
import sqlite3
import conf
import os


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


def new_db_connection(path2db):
    """
    :param path2db: (string) path where the database will be or is located
    """
    return sqlite3.connect(path2db + "fasterDB_lite.db")


def creation_of_gene_table(new_db):
    """
    Create a genes table in ``new_db``

    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE genes (
        id int(10),
        official_symbol VARCHAR(17) NOT NULL,
        chromosome VARCHAR(2) NOT NULL,
        strand tinyint(1) NOT NULL,
        start_on_chromosome int(10) NOT NULL,
        end_on_chromosome int(10) NOT NULL,
        sequence mediumtext NOT NULL,
        PRIMARY KEY (id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def creation_of_intron_table(new_db):
    """
    Create an intron_genomiques table in ``new_db``

    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE intron_genomiques(
        id_gene int(10),
        pos_sur_gene int(10),
        start_on_gene int(10) NOT NULL,
        end_on_gene int(10) NOT NULL,
        PRIMARY KEY(id_gene, pos_sur_gene),
        FOREIGN KEY (id_gene) REFERENCES genes(id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def creation_of_exon_table(new_db):
    """
    Create an exon_genomiques table in ``new_db``

    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE exon_genomiques (
        id_gene int(10) NOT NULL,
        pos_sur_gene int(10) NOT NULL,
        start_on_gene int(10) NOT NULL,
        end_on_gene int(10) NOT NULL,
        exon_type VARCHAR(3),
        cds_start_on_gene int(10),
        cds_end_on_gene int(10),
        offset_before_exon tinyint(2),
        offset_after_exon tinyint(2),
        PRIMARY KEY(id_gene, pos_sur_gene),
        FOREIGN KEY (id_gene) REFERENCES genes(id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def database_creator():
    """
    Create an empty database
    """
    out_path = os.path.realpath(__file__)
    out_path = "/".join(out_path.split("/")[:-2]) + "/result/"
    if not os.path.isdir(out_path):
        out_path = os.path.dirname(os.path.realpath(__file__))
    new_db = new_db_connection(out_path)
    creation_of_gene_table(new_db)
    creation_of_intron_table(new_db)
    creation_of_exon_table(new_db)


if __name__ == "__main__":
    database_creator()
