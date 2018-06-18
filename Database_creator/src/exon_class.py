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
        self.gene = Gene(gene_name, gene_id)
        self.position = exon_position
        length, exon_type, donor, acceptor = self.get_exon_info(cnx)
        self.length = length
        self.acceptor = acceptor
        self.donor = donor
        self.type = exon_type

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
        WHERE id_gene = \"""" + str(self.gene.id) + """\"
        AND pos_on_gene = """ + str(self.position) + """ ;
        """
        cursor.execute(query)
        result = cursor.fetchone()
        if result is None:
            return None, None, None, None
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
        printd("Exon " + str(gene_name) + "_" + str(exon_position))
        ExonClassMain.__init__(self, cnx, gene_name, gene_id, exon_position)
        self.gene.gene_filler(cnx)
        self.upstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position - 1)
        self.downstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position + 1)
        self.upstream_intron = Intron(cnx, gene_id, self.position - 1, self.gene.sequence, "upstream")
        self.downstream_intron = Intron(cnx, gene_id, self.position, self.gene.sequence, "downstream")
        iupac, dnt = self.get_iupac_dnt_exon(cnx)
        self.iupac = iupac
        self.dnt = dnt
        # once the exon is fully created we delete the gene sequence for memory efficiency
        self.gene.sequence = None

    def get_iupac_dnt_exon(self, cnx):
        """
        Get the iupac and dnt content of the exon
        :return: (string) the iupac content of the exon
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
            return None, None
        sequence = self.gene.sequence[start:stop]
        printd("Exon sequence:")
        printd(sequence)
        iupac = iupac_frequencies(sequence)
        dnt = dinucleotide_frequencies(sequence)
        return ";".join(list(map(str, iupac))), ";".join(list(map(str, dnt)))


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

    def get_iupac_dnt_and_gene_length(self, cnx):
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
        self.length = len(sequence)
        iupac = iupac_frequencies(sequence)
        res = ";".join(list(map(str, iupac)))
        self.iupac = res
        dnt = dinucleotide_frequencies(sequence)
        res = ";".join(list(map(str, dnt)))
        self.dnt = res

    def get_nb_intron_and_median_intron_size(self, cnx):
        """
        Get the number of intron in a gene and the median intron size

        :param cnx: (sqlite3 object) connection to fasterDB
        """
        cursor = cnx.cursor()
        query = """SELECT end_on_gene - start_on_gene + 1
                    FROM introns
                    WHERE id_gene = """ + str(self.id) + ";"
        cursor.execute(query)
        result = cursor.fetchall()
        self.nb_intron = len(result)
        if not result:
            self.median_intron_size = None
        else:
            self.median_intron_size = int(np.median([size[0] for size in result]))

    def gene_filler(self, cnx):
        """
        Executes the functions ``get_iupac_and_gene_length`` and \
        ``get_nb_intron_and_median_intron_size``
        :param cnx:(sqlite3 object) connection to fasterDB
        :return:
        """
        self.get_iupac_dnt_and_gene_length(cnx)
        self.get_nb_intron_and_median_intron_size(cnx)


class Intron:
    """Create an intron"""
    def __init__(self, cnx, gene_id, pos_on_gene, gene_sequence, location):
        """

        :param cnx: (sqlite3 object) allow connection to **FasterDB Lite** database
        :param gene_id:  (int) the id of a gene
        :param pos_on_gene: (int) an intron position on a gene
        :param gene_sequence: (string) the sequence of the gene, gene_id
        :param location: (string) upstream or downstream
        """
        self.gene_id = gene_id
        self.position = pos_on_gene
        self.location = location
        self.length = None
        self.iupac_proxi = None
        self.iupac_dist = None
        self.dnt_proxi = None
        self.dnt_dist = None
        self.get_intron_information(cnx, gene_sequence)

    def get_intron_information(self, cnx, gene_seq):
        """
        Get the intron size and the iupac composition of distal and proximal sequence

        :param cnx: (sqlite3 object) allows connection to **FasterDB Lite** database
        :param gene_seq: (string) the sequence of the gene, gene_id
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
            return None, None, None
        start = result[0][0] - 1
        end = result[0][1]
        sequence = gene_seq[start:end]
        printd("Intron " + self.location + " sequence : ")
        printd(sequence)
        if self.location == "upstream":
            sequence = sequence[::-1]
        proxi_seq = sequence[0:25]
        distal_seq = sequence[26:101]
        if len(proxi_seq) > 0:
            self.iupac_proxi = ";".join(list(map(str, iupac_frequencies(proxi_seq))))
            self.dnt_proxi = ";".join(list(map(str, dinucleotide_frequencies(proxi_seq))))
        else:
            self.iupac_proxi = None
            self.dnt_proxi = None
        if len(distal_seq) > 0:
            self.iupac_dist = ";".join(list(map(str, iupac_frequencies(distal_seq))))
            self.dnt_dist = ";".join(list(map(str, dinucleotide_frequencies(distal_seq))))
        else:
            self.iupac_dist = None
            self.dnt_dist = None
        self.length = len(sequence)


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


def dinucleotide_frequency(sequence, dnt):
    """
    Get the frequency (in %) of ``dnt`` in ``sequence``

    :param sequence: (string) a nucleotide sequence
    :param dnt: (string) a di-nucleotide
    :return: (float) the frequency of the di-nucleotide ``dnt`` in ``sequence``
    """
    seqlen = len(sequence) - (len(dnt) - 1)
    count = 0
    dnt1 = dnt[0]
    dnt2 = dnt[1]
    for i in range(seqlen):
        if sequence[i] == dnt1 and sequence[i + 1] == dnt2:
            count += 1
    return round(float(count) / seqlen * 100, 1)


def dinucleotide_frequencies(sequence):
    """
    Get di-nucleotides frequencies for a sequence.

    :param sequence: (string) a nucleotide sequence
    :return:  (list of float) the frequency of di-nucleotides \
    AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively
    """
    dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    results = []
    for dnt in dnt_list:
        results.append(dinucleotide_frequency(sequence, dnt))
    return results


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
