#!/usr/bin/env python3

# -*- coding: uft-8 -*-

"""
Description:

    This script allows to retrieve every nucleotide (iupac) frequencies within \
    the exons, and 100 bp before and after the exons.
"""


class ExonClass:
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
        printd("Exon " + str(gene_name) + "_" + str(exon_position))
        self.gene = Gene(cnx, gene_name, gene_id)
        self.position = exon_position
        self.upstream_intron = None
        self.downstream_intron = None
        self.exon = None
        self.get_exon_data(cnx)

        # calculating iupac frequencies of exons, upstream intron and downstream_intron
        self.iupac_exon = iupac_frequencies(self.exon)
        self.iupac_up_intron = iupac_frequencies(self.upstream_intron)
        self.iupac_down_intron = iupac_frequencies(self.downstream_intron)

        # debug display
        printd("upstream intron sequence : ")
        printd(self.upstream_intron)
        printd("upstream intron iupac : ")
        printd(self.iupac_up_intron)
        printd("exon sequence : ")
        printd(self.exon)
        printd("exon iupac : ")
        printd(self.iupac_exon)
        printd("downstream intron sequence : ")
        printd(self.downstream_intron)
        printd("downstream intron iupac : ")
        printd(self.iupac_down_intron)

        del self.upstream_intron
        del self.exon
        del self.downstream_intron
        del self.gene.sequence

    def get_exon_data(self, cnx):
        """
        Get the sequences of exons, upstream and downstream intron (100bp before and after exon).

        :param cnx: (sqlite3 object) connection to fasterDB
        """
        cursor = cnx.cursor()
        query = """SELECT start_on_gene, end_on_gene
                   FROM exons
                   WHERE id_gene = """ + str(self.gene.id) + """
                   AND pos_on_gene = """ + str(self.position) + ";"
        cursor.execute(query)
        result = cursor.fetchall()
        if len(result) > 1:
            print("More than one exon was retrieve with a sing exon id")
            print("Exiting...")
            exit(1)
        start = result[0][0] - 1
        stop = result[0][1]
        if start >= stop:
            self.exon = None
        self.exon = self.gene.sequence[start:stop]
        if self.position != 1:
            query = """SELECT end_on_gene
                       FROM exons
                       WHERE id_gene = """ + str(self.gene.id) + """
                       AND pos_on_gene = """ + str(self.position - 1) + ";"
            cursor.execute(query)
            result = cursor.fetchall()
            end_previous = result[0][0] + 1
            if end_previous-1 >= start:
                self.upstream_intron = None
            else:
                self.upstream_intron = self.gene.sequence[end_previous-1:start]
        else:
            self.upstream_intron = None

        query = """SELECT start_on_gene
                   FROM exons
                   WHERE id_gene = """ + str(self.gene.id) + """
                   AND pos_on_gene = """ + str(self.position + 1) + ";"
        cursor.execute(query)
        result = cursor.fetchall()
        if result:
            start_next = result[0][0] - 2
            if stop >= start_next+1:
                self.downstream_intron = None
            else:
                self.downstream_intron = self.gene.sequence[stop:start_next+1]
        else:
            self.downstream_intron = None


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
        self.sequence = None
        self.get_sequence_gene(cnx)

    def get_sequence_gene(self, cnx):
        """
        Get gene sequence

        :param cnx: (sqlite3 object) connection to fasterDB
        """
        cursor = cnx.cursor()
        query = """SELECT sequence
                   FROM genes
                   WHERE id = """ + str(self.id) + """;"""
        cursor.execute(query)
        result = cursor.fetchall()
        if len(result) > 1:
            print("Multiple sequence were retrieve from a single gene id")
            print("Exiting...")
            exit(1)
        sequence = result[0][0]
        self.sequence = sequence


def iupac_frequencies(sequence):
    """
    Get iupac-nucleotides frequencies for a sequence.

    :param sequence: (string) a hexa-nucleotide sequence
    :return:  (list of float) the frequency of hexa-nucleotides
    """
    iupac_dic = {"R": ["A", "G"], "Y": ["C", "T"], "W": ["A", "T"], "S": ["C", "G"], "K": ["T", "G"], "M": ["C", "A"]}
    res_dic = {}
    if sequence is None:
        return None
    for nt in ["A", "C", "G", "T"]:
        res_dic[nt] = float(sequence.count(nt)) / len(sequence) * 100
    for nt in iupac_dic.keys():
        res_dic[nt] = res_dic[iupac_dic[nt][0]] + res_dic[iupac_dic[nt][1]]
    return res_dic


def set_debug(debug=0):
    """
    Set debug mod if ``debug`` == 1

    :param debug: (int) 0 no debug, 1 debug mode
    """
    global d
    d = debug


def printd(message):
    """

    :param message: (string) message we want to print
    :return:
    """
    if d == 1:
        print(message)
