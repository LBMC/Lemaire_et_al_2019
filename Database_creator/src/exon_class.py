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
        else:
            if result[0] > 0:
                length_exon = result[0]
            else:
                length_exon = None
        return length_exon, result[1], result[2], result[3]


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
        iupac, dnt, exon_sequence = self.get_iupac_dnt_exon(cnx)
        self.iupac = iupac
        self.dnt = dnt
        self.upstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position - 1)
        self.downstream_exon = ExonClassMain(cnx, gene_name, gene_id, exon_position + 1)
        self.upstream_intron = Intron(cnx, gene_id, self.position - 1, self.gene.sequence, "upstream", exon_sequence)
        self.downstream_intron = Intron(cnx, gene_id, self.position, self.gene.sequence, "downstream", exon_sequence)
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
        if not full_defined(sequence):
            return None, None
        printd("Exon sequence:")
        printd(sequence)
        iupac = iupac_frequencies(sequence)
        if len(sequence) > 1:
            dnt = ";".join(list(map(str, dinucleotide_frequencies(sequence))))
        else:
            dnt = None
        return ";".join(list(map(str, iupac))), dnt, sequence


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
        if len(sequence) > 0:
            self.length = len(sequence)
        else:
            self.length = None
        if full_defined(sequence):
            iupac = iupac_frequencies(sequence)
            res = ";".join(list(map(str, iupac)))
            self.iupac = res
            dnt = dinucleotide_frequencies(sequence)
            res = ";".join(list(map(str, dnt)))
            self.dnt = res
        else:
            self.iupac = None
            self.dnt = None

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
            self.median_intron_size = int(np.median([size[0] for size in result if size[0] > 0]))

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
    def __init__(self, cnx, gene_id, pos_on_gene, gene_sequence, location, exon_sequence):
        """

        :param cnx: (sqlite3 object) allow connection to **FasterDB Lite** database
        :param gene_id:  (int) the id of a gene
        :param pos_on_gene: (int) an intron position on a gene
        :param gene_sequence: (string) the sequence of the gene, gene_id
        :param location: (string) upstream or downstream
        :param exon_sequence: (string) the sequence of the exon of interest
        """
        self.gene_id = gene_id
        self.position = pos_on_gene
        self.location = location
        self.length = None
        self.iupac_proxi = None
        self.iupac = None
        self.iupac_ei = None  # Iupac 100 nt intron 50 nt exons
        self.dnt_proxi = None
        self.dnt = None
        self.dnt_ei = None
        self.get_intron_information(cnx, gene_sequence, exon_sequence)

    def get_intron_information(self, cnx, gene_seq, exon_sequence):
        """
        Get the intron size and the iupac composition of distal and proximal sequence

        :param cnx: (sqlite3 object) allows connection to **FasterDB Lite** database
        :param gene_seq: (string) the sequence of the gene, gene_id
        :param exon_sequence: (string) the sequence of the interest exons
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
        if len(sequence) > 0 and full_defined(sequence):
            self.iupac = ";".join(list(map(str, iupac_frequencies(sequence))))
        if len(sequence) > 1 and full_defined(sequence):
            self.dnt = ";".join(list(map(str, dinucleotide_frequencies(sequence))))
        if self.location == "upstream":
            sequence = sequence[::-1]
        proxi_seq = sequence[0:100]
        if self.location == "upstream":
            proxi_seq = proxi_seq[::-1]
        if len(proxi_seq) > 0 and full_defined(proxi_seq):
            self.iupac_proxi = ";".join(list(map(str, iupac_frequencies(proxi_seq)))) + ";%s" % (len(proxi_seq))
        if len(proxi_seq) > 1 and full_defined(proxi_seq):
            self.dnt_proxi = ";".join(list(map(str, dinucleotide_frequencies(proxi_seq)))) + ";%s" % (len(proxi_seq))
        if len(proxi_seq) > 0 and len(exon_sequence) > 0:
            if self.location == "upstream":
                exon_intron = proxi_seq + exon_sequence[0:50]
            else:
                exon_intron = exon_sequence[-50:] + proxi_seq
            if full_defined(exon_intron):
                self.iupac_ei = ";".join(list(map(str, iupac_frequencies(exon_intron)))) + ";%s" % (len(exon_intron))
                self.dnt_ei = ";".join(list(map(str, dinucleotide_frequencies(exon_intron)))) + \
                              ";%s" % (len(exon_intron))
            printd("Exon intron %s jonction sequence" % self.location)
            printd(exon_intron)
        if len(sequence) > 0:
            self.length = len(sequence)
        else:
            self.length = None


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
    seq_len = sequence.count("A") + sequence.count("T") + sequence.count("G") + sequence.count("C")
    for nt in ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]:
        if nt not in iupac:
            result.append(round((float(sequence.count(nt)) / seq_len) * 100, 1))
        else:
            result.append(round((float(sequence.count(iupac[nt][0]) +
                                       sequence.count(iupac[nt][1])) / seq_len) * 100, 1))
    return result


def dinucleotide_frequencies(sequence):
    """
    Get di-nucleotides frequencies for a sequence.

    :param sequence: (string) a nucleotide sequence
    :return:  (list of float) the frequency of di-nucleotides \
    AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively
    """

    dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    dnt_dic = {"AA": 0., "AC": 0., "AG": 0., "AT": 0., "CA": 0., "CC": 0., "CG": 0., "CT": 0., "GA": 0.,
               "GC": 0., "GG": 0., "GT": 0., "TA": 0., "TC": 0., "TG": 0., "TT": 0.}
    results = []
    seqlen = len(sequence) - 1
    count = 0
    for i in range(seqlen):
        cdnt = sequence[i] + sequence[i+1]
        if cdnt in dnt_dic:
            dnt_dic[cdnt] += 1
            count += 1
    for dnt in dnt_list:
        results.append(round(dnt_dic[dnt] / count * 100, 1))
    return results


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
