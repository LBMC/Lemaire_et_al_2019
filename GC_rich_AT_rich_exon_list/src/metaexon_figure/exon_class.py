#!/usr/bin/python3.5

"""
Description:

    This script contains an exon class that will allow to get, \
    For each exons their nucleotides frequencies in the regions [-100;+50] and [-50+100]  \
    surrounding the exons (frequencies calculated using a spliding windows
    database.

"""


import math
import numpy as np


class ExonClassMain:
    """
    A class corresponding to an exon. This class allows to \
    extract easily information about an exons.
    """
    def __init__(self, cnx, gene_name, gene_id, exon_position):
        """
        Initiate the creation of an exon

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of a gene
        :param exon_position:  (int) the position of the exon on the gene
        """
        self.gene = Gene(gene_name, gene_id)
        self.position = exon_position
        length, exon_type, start, end = self.get_exon_info(cnx)
        self.length = length
        self.type = exon_type
        self.start = start
        self.end = end

    def get_exon_info(self, cnx):
        """
        Get basics data about an exon.

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :return: (int) the length of the exon
        """
        cursor = cnx.cursor()
        query = """
        SELECT end_on_gene - start_on_gene + 1, exon_type, start_on_gene, end_on_gene
        FROM exons
        WHERE id_gene = \"""" + str(self.gene.id) + """\"
        AND pos_on_gene = """ + str(self.position) + """ ;
        """
        cursor.execute(query)
        result = cursor.fetchone()
        if result is None:
            return None, None
        else:
            if result[0] > 0:
                length_exon = result[0]
            else:
                length_exon = None
        return length_exon, result[1], result[2] - 1, result[3]


class ExonClass(ExonClassMain):
    """
    Contains every data of interest
    """
    def __init__(self, cnx, gene_name, gene_id, exon_position, window_size):
        """
        Initiate the creation of an exon

        :param cnx: (sqlite3 object) connection to fasterDB Lite
        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of a gene
        :param exon_position:  (int) the position of the exon on the gen
        :param window_size: (size of the window for the metaexon
        """
        printd("Exon " + str(gene_name) + "_" + str(exon_position))
        ExonClassMain.__init__(self, cnx, gene_name, gene_id, exon_position)
        self.gene.gene_filler(cnx)
        self.upstream_intron = Intron(cnx, gene_id, self.position - 1, "upstream")
        self.downstream_intron = Intron(cnx, gene_id, self.position, "downstream")
        self.sequence_5p_ext = None
        self.sequence_3p_ext = None
        self.retrieve_adjacent_exon_sequence(window_size)
        printd("sequence_5p_ext")
        printd(self.sequence_5p_ext)
        printd("sequence_3p_ext")
        printd(self.sequence_3p_ext)
        # once the exon is fully created we delete the gene sequence for memory efficiency
        self.gene.sequence = None

    def retrieve_adjacent_exon_sequence(self, window_size):
        """
        Retrieve the adjacent exon sequences
        :param window_size: (int) the size of the windows
        :return: (2 string) the sequences surrounding the introns
        """
        val_start = 100 + math.floor(window_size/2)
        val_stop = 50 + round(window_size/2)
        if self.upstream_intron.start is not None and self.upstream_intron.end \
                and self.start is not None and self.end is not None:
            upstream_sequence = self.gene.sequence[self.upstream_intron.start:self.upstream_intron.end]
            exon_sequence = self.gene.sequence[self.start:self.end]
            if len(upstream_sequence) > val_start:
                upstream_sequence = upstream_sequence[-val_start:]
            else:
                nb_n = abs(len(upstream_sequence) - val_start)
                upstream_sequence = "N" * nb_n + upstream_sequence
            if len(exon_sequence) > val_stop:
                exon_sequence = exon_sequence[0:val_stop]
            else:
                nb_n = abs(len(exon_sequence) - val_stop)
                exon_sequence += "N" * nb_n
            sequence_5p_ext = upstream_sequence + exon_sequence

        else:
            sequence_5p_ext = None

        if self.downstream_intron.start is not None and self.downstream_intron.end and \
                self.start is not None and self.end is not None:
            downstream_sequence = self.gene.sequence[self.downstream_intron.start:self.downstream_intron.end]
            exon_sequence = self.gene.sequence[self.start:self.end]
            if len(exon_sequence) > val_stop:
                exon_sequence = exon_sequence[-val_stop:]
            else:
                nb_n = abs(len(exon_sequence) - val_stop)
                exon_sequence = "N" * nb_n + exon_sequence

            if len(downstream_sequence) > val_start:
                downstream_sequence = downstream_sequence[0:val_start]
            else:
                nb_n = abs(len(downstream_sequence) - val_start)
                downstream_sequence += "N" * nb_n
            sequence_3p_ext = exon_sequence + downstream_sequence
        else:
            sequence_3p_ext = None

        self.sequence_5p_ext = sequence_5p_ext
        self.sequence_3p_ext = sequence_3p_ext


class Gene:
    """
    Contains every info of interest for a gene
    """
    def __init__(self, gene_name, gene_id):
        """
        Initiate the creation of a gene

        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of the gene
        """
        self.name = gene_name
        self.id = gene_id
        self.nb_intron = None
        self.median_intron_size = None
        self.iupac = None
        self.dnt = None
        self.length = None
        self.sequence = None

    def gene_filler(self, cnx):
        """
        Get iupac-dnt information and length of a gene

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
        if len(sequence) > 0:
            self.length = len(sequence)
        else:
            self.length = None


class Intron:
    """Create an intron"""
    def __init__(self, cnx, gene_id, pos_on_gene, location):
        """

        :param cnx: (sqlite3 object) allow connection to **FasterDB Lite** database
        :param gene_id:  (int) the id of a gene
        :param pos_on_gene: (int) an intron position on a gene
        :param location: (string) upstream or downstream
        """
        self.gene_id = gene_id
        self.position = pos_on_gene
        self.location = location
        self.start = None
        self.end = None
        self.get_intron_information(cnx)

    def get_intron_information(self, cnx):
        """
        Get the intron coordinates on it's gene

        :param cnx: (sqlite3 object) allows connection to **FasterDB Lite** database
        :return: The following information is returned
            - length: (int) intron length
            - proximal: (string) the proximal (0:25) frequency in the intron
            - distal: (string)  the proximal (26:100) iupac frequency in the intorn
        """
        cursor = cnx.cursor()
        query = """SELECT start_on_gene, end_on_gene
                   FROM introns
                   WHERE id_gene = """ + str(self.gene_id) + """
                   AND pos_on_gene = """ + str(self.position) + ";"
        cursor.execute(query)
        result = cursor.fetchall()
        if len(result) > 1:
            print("Error : multiple intron find from a simple intron id")
            print("exiting...")
            exit(1)
        if not result:
            return None, None
        start = result[0][0] - 1
        end = result[0][1]
        self.start = start
        self.end = end


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


def full_defined(sequence):
    """
    Says if all the nucleotide are well defined within the sequence

    :param sequence: (string) nucleotide sequence
    :return: (boolean) True if all nucleotide are well defined, False else
    """
    seq_defined = sequence.count("A") + sequence.count("T") + sequence.count("G") + sequence.count("C")
    if seq_defined / len(sequence) >= 0.95:
        return True
    return False


def get_metagene_vectors_windowsed(exon_list, window_size):
    """

    :param exon_list: (list of list of 2 ints) list of exons identified by their gene_id and exon position
    :param window_size: (int) the size of the window we have choose for the metagene analysis \
    :return: 4 variables:

    * final_res_5p : (dict of list of 150 floats) for each amino acid (keys of the dictionary) gives its mean \
    proportion in a windows of 'window_size' size in the sequence 5p \
    (100 nt before the exon and the 50 first nt of this exon) for every exon of interest. The the value of \
    this window is averaged with the values of every studied exons (same windows for the same sequence at the same \
    position of this sequence).
    * final_res_3 : (dict of list of 150 floats) for each amino acid (keys of the dictionary) gives its mean \
    proportion in a windows of 'window_size' size in the sequence 3p \
    (50 last nt of a studied exon and the 100 nt after this exon) \
    for every exon of interest. The the value of \
    this window is averaged with the values of every studied exons (same windows for the same sequence at the same \
    position of this sequence).
    * p5_analyzed (int) the number of sequences 5p analyzed
    * p3_analyzed (int) the number of sequences 3p analyzed
    """
    p5_analyzed = 0
    p3_analyzed = 0
    i = 0
    if exon_list[i].sequence_5p_ext is not None:
        res_5p = [{'A': [], 'C': [], 'G': [], 'T': []} for i in range(len(exon_list[i].sequence_5p_ext) - window_size)]
        res_3p = [{'A': [], 'C': [], 'G': [], 'T': []} for i in range(len(exon_list[i].sequence_5p_ext) - window_size)]
    while exon_list[i].sequence_5p_ext is None:
        res_5p = [{'A': [], 'C': [], 'G': [], 'T': []} for i in range(len(exon_list[i].sequence_5p_ext) - window_size)]
        res_3p = [{'A': [], 'C': [], 'G': [], 'T': []} for i in range(len(exon_list[i].sequence_5p_ext) - window_size)]
        i += 1

    for i in range(len(exon_list)):
        if exon_list[i].sequence_5p_ext is not None:
            p5_analyzed += 1
            for j in range(len(exon_list[i].sequence_5p_ext) - window_size):
                for key in res_5p[j].keys():
                    cur_seq = exon_list[i].sequence_5p_ext[j:j + window_size]
                    if full_defined(cur_seq):
                        seq_len = cur_seq.count("A") + cur_seq.count("C") + cur_seq.count("G") + cur_seq.count("T")
                        res_5p[j][key].append(float(cur_seq.count(key)) / seq_len)
                    else:
                        res_5p[j][key].append(None)
        if exon_list[i].sequence_3p_ext is not None:
            p3_analyzed += 1
            for j in range(len(exon_list[i].sequence_3p_ext) - window_size):
                for key in res_3p[j].keys():
                    cur_seq = exon_list[i].sequence_3p_ext[j:j + window_size]
                    if full_defined(cur_seq):
                        seq_len = cur_seq.count("A") + cur_seq.count("C") + cur_seq.count("G") + cur_seq.count("T")
                        res_3p[j][key].append(float(cur_seq.count(key)) / seq_len)

    final_res_5p = {'A': [], 'C': [], 'G': [], 'T': []}
    final_res_3p = {'A': [], 'C': [], 'G': [], 'T': []}
    for i in range(len(res_5p)):
        for key in final_res_5p.keys():
            final_res_5p[key].append(round(np.nanmean(np.array(res_5p[i][key], dtype=float)) * 100, 4))
            final_res_3p[key].append(round(np.nanmean(np.array(res_3p[i][key], dtype=float)) * 100, 4))
    return final_res_5p, final_res_3p, p5_analyzed, p3_analyzed
