#!/usr/bin/python3.5

"""
Description:

    This script contains an exon class that will allow to get, \
    For each exons numbers of interest information from **FasterDB Lite** \
    database.

"""

import numpy as np


class ExonClassMain:
    """
    A class corresponding to an exon. This class allows to \
    extract easily information about and exons and its vicinity sequence.
    """
    def __init__(self, cnx, gene_name, gene_id, exon_position):
        """
        Initiate the creation of an exon

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of a gene
        :param exon_position:  (int) the position of the exon on the gene
        """
        self.gene = Gene(cnx, gene_name, gene_id)
        self.position = exon_position
        length, exon_type, donor, acceptor = self.get_exon_info(cnx)
        self.length = length
        self.acceptor = acceptor
        self.donor = donor
        self.exon_type = exon_type

    def get_exon_info(self, cnx):
        """
        Check if an exon exist

        :param cnx: (sqlite3 object) connection to fasterDB Lite
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
        result = cursor.fetchone()
        if result is None:
            return (None, None, None, None)
        return result


class ExonClass(ExonClassMain):
    """
    Contains every data of interest
    """
    def __init__(self, cnx, gene_name, gene_id, exon_position):
        """
        Initiate the creation of an exon

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of a gene
        :param exon_position:  (int) the position of the exon on the gene
        """
        ExonClassMain.__init__(self, cnx, gene_name, gene_id, exon_position)
        self.upstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position - 1)
        self.downstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position + 1)


class Gene:
    """
    Contains every info of interest for a gene
    """
    def __init__(self, cnx, gene_name, gene_id):
        """
        Initiate the creation of a gene

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of the gene
        """
        self.name = gene_name
        self.id = gene_id
        self.nb_intron = None
        self.median_intron_size = None
        self.iupac = None
        self.length = None

    def get_iupac_and_gene_length(self, cnx):
        """
        Get iupac information and length of a gene

        :param cnx: (sqlite3 object) connection to fasterDB
        """
        cursor = cnx.cursor()
        query = """SELECT sequence
                   FROM genes
                   WHERE id = """ + self.id + """;"""
        cursor.execute(query)
        result = cursor.fetchall()
        if len(result) > 1:
            print("Multiple sequence were retrive from a single gene id")
            print("Exiting...")
            exit(1)
        sequence = result[0][0]
        self.length = len(sequence)
        iupac = iupac_frequencies(sequence)
        res = ";".join(list(map(str,iupac)))
        self.iupac = res

    def get_nb_intron_and_median_intron_size(self, cnx):
        """
        Get the number of intron in a gene and the median intron size

        :param cnx: (sqlite3 object)
        """
        cursor = cnx.cursor()
        query = """SELECT end_on_gene - start_on_gene + 1
                    FROM introns
                    WHERE id_gene = """ + self.id + ";"
        cursor.execute(query)
        result = cursor.fetchall()
        self.nb_intron = len(result)
        self.median_intron_size = np.median([size[0] for size in result])


# simple function for getting iupac frequencies
def iupac_frequencies(sequence):
    """
    Get iupac frequencies info for a sequence

    :param sequence: (string) a nucleotide sequence
    :return: (list of float) the frequency of nucleotides A, C, G, T, S, W, R, Y, K, M respectively
    """
    iupac = {"S": ["C", "G"], "W": ["A", "T"], "R": ["A", "G"],
             "Y": ["C", "T"], "K": ["T", "G"], "M": ["A", "C"]}
    result = []
    for nt in ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]:
        if nt not in iupac:
            result.append(round((float(sequence.count(nt)) / len(sequence)) * 100, 1))
        else:
            result.append(round((float(sequence.count(iupac[nt][0]) +
                                       sequence.count(iupac[nt][1])) / len(sequence)) * 100, 1))
    return result

