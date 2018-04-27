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
        length, exon_type, donor, acceptor = self.get_exon_info(cnx)
        self.length = length
        self.acceptor = acceptor
        self.donor = donor
        self.exon_type = exon_type

    def get_exon_info(self, cnx):
        """
        Check if an exon exist

        :param cnx: (pymysql object) connection to fasterDB
        :return: (int) the length of the exon
        """
        cursor = cnx.cursor()
        query = """
        SELECT end_on_gene - start_on_gene + 1, exon_type, force_donor, force_acceptor
        FROM exons
        WHERE id_gene = \"""" + str(self.gene_id) + """\"
        AND pos_on_gene = """ + str(self.position) + """ ;
        """
        cursor.execute(query)
        if cursor.arraysize == 0:
            return (None, None, None, None)
        return cursor.fetchone()


