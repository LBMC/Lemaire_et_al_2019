#!/usr/bin/python3.5

"""
Description:

    This script contains an exon class that will allow to get, \
    For each exons numbers of interest information from **FasterDB Lite** \
    database.

"""


class ExonClassMain:
    """
    A class corresponding to an exon. This class allows to \
    extract easily information about and exons and its vicinity sequence.
    """
    def __init__(self, cnx, gene_name, gene_id, exon_position):
        """
        Initiate the creation of an exon

        :param cnx: (pymysql object) connection to fasterDB
        :param gene_name: (string) official symbol for an exons
        :param exon_position:  (int) the position of the exon on the gene
        """
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.position = exon_position
        self.length = self.get_exon_length(cnx)
        self.acceptor = self.get_force(cnx, 0)
        self.donor = self.get_force(cnx, 1)

    def get_exon_length(self, cnx):
        """
        Check if an exon exist

        :param cnx: (pymysql object) connection to fasterDB
        :return: (int) the length of the exon
        """
        cursor = cnx.cursor()
        query = """
        SELECT end_on_gene - start_on_gene + 1
        FROM exons
        WHERE id_gene = \"""" + str(self.gene_id) + """\"
        AND pos_on_gene = """ + str(self.position) + """ ;
        """
        cde = cursor.execute(query)
        if cde == 0:
            return None
        return cursor.fetchone()[0]

    def get_force(self, cnx, isdonor):
        """
        Get the force of the donor or the acceptor of an exons

        :param cnx: (pymysql object) connection to fasterDB
        :param isdonor: (int) 1 for donor 0 for acceptor
        """
        cursor = cnx.cursor()
        if isdonor == 1:
            column = "force_donor"
        else:
            column = "force_acceptor"
        query = """
        SELECT """ + column + """
        FROM exons
        WHERE id_gene = """ + str(self.gene_id) + """
        AND pos_on_gene = """ + str(self.position) + """ ; """
        cde = cursor.execute(query)
        if cde == 0:
            return None
        result = cursor.fetchall()[0]
        return result
