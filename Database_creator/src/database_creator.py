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
import sqlite3
import conf
import os

base_name = "fasterDB_lite.db"
out_path = "/".join(os.path.realpath(__file__).split("/")[:-2]) +  "/result/"


# Functions
def new_db_connection(file_name):
    """
    :param file_name: (string) file where the databse will be created
    """
    return sqlite3.connect(file_name)


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
    CREATE TABLE introns(
        id_gene int(10),
        pos_on_gene int(10),
        start_on_gene int(10) NOT NULL,
        end_on_gene int(10) NOT NULL,
        start_on_chromosome int(10) NOT NULL,
        end_on_chromosome int(10) NOT NULL,
        PRIMARY KEY(id_gene, pos_on_gene),
        FOREIGN KEY (id_gene) REFERENCES genes(id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def creation_of_exon_table(new_db):
    """
    Create an exon_partial table in ``new_db``

    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE exon_partial (
        id_gene int(10) NOT NULL,
        pos_on_gene int(10) NOT NULL,
        start_on_gene int(10) NOT NULL,
        end_on_gene int(10) NOT NULL,
        exon_type VARCHAR(3),
        cds_start_on_gene int(10),
        cds_end_on_gene int(10),
        offset_before_exon tinyint(2),
        offset_after_exon tinyint(2),
        chromosome VARCHAR(2),
        start_on_chromosome INT(12),
        end_on_chromosome INT(12),
        PRIMARY KEY(id_gene, pos_on_gene),
        FOREIGN KEY (id_gene) REFERENCES genes(id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def creation_of_full_exon_table(new_db):
    """
    Create an exon_genomiques table in ``new_db``

    :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
    """
    cursor = new_db.cursor()
    query = """
    CREATE TABLE exons (
        id_gene int(10) NOT NULL,
        pos_on_gene int(10) NOT NULL,
        start_on_gene int(10) NOT NULL,
        end_on_gene int(10) NOT NULL,
        exon_type VARCHAR(3),
        cds_start_on_gene int(10),
        cds_end_on_gene int(10),
        offset_before_exon tinyint(2),
        offset_after_exon tinyint(2),
        force_donor int,
        force_acceptor int,
        chromosome VARCHAR(2),
        start_on_chromosome INT(12),
        end_on_chromosome INT(12),
        PRIMARY KEY(id_gene, pos_on_gene),
        FOREIGN KEY (id_gene) REFERENCES genes(id)
    );
    """
    cursor.execute(query)
    new_db.commit()


def creation_of_force_splicing_table(new_db):
    """
       Create a force splicing site table in ``new_db``

       :param new_db: (sqlite3 object) all the info we need to connect to sqlite3
       """
    cursor = new_db.cursor()
    query = """
       CREATE TABLE force_splicing_site (
           id_gene int(10) NOT NULL,
           is_donor tinyint(1),
           pos_on_gene int(10) NOT NULL,
           force float,
           FOREIGN KEY(id_gene, pos_on_gene) REFERENCES exon_partial(id_gene, pos_on_gene),
           FOREIGN KEY (id_gene) REFERENCES genes(id)
       );
       """
    cursor.execute(query)
    new_db.commit()


def database_creator():
    """
    Create an empty database
    """
    new_db = new_db_connection(out_path + base_name)
    creation_of_gene_table(new_db)
    creation_of_intron_table(new_db)
    creation_of_exon_table(new_db)
    creation_of_full_exon_table(new_db)
    creation_of_force_splicing_table(new_db)
    new_db.close()
    return out_path + base_name


if __name__ == "__main__":
    database_creator()
