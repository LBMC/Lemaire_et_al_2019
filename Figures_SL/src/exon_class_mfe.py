#!/usr/bin/python3.5

"""
Description:

    This script contains an exon class that will allow to get, \
    For each exons their upstream intronic sequence

"""


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
        self.start = None
        self.stop = None
        self.get_exon_position(cnx)

    def get_exon_position(self, cnx):
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
        self.start = result[0][0] - 1
        self.stop = result[0][1]


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
        :param exon_position:  (int) the position of the exon on the gen
        """
        printd("Exon " + str(gene_name) + "_" + str(exon_position))
        ExonClassMain.__init__(self, cnx, gene_name, gene_id, exon_position)
        self.gene.gene_filler(cnx)
        self.seq_3ss = None
        self.seq_5ss = None
        if self.start - 25 >= 0 and self.start + 25 < len(self.gene.sequence):
            self.seq_3ss = self.gene.sequence[self.start - 25: self.start + 25]
        if self.stop - 25 >= 0 and self.stop + 25 < len(self.gene.sequence):
            self.seq_5ss = self.gene.sequence[self.stop -25: self.stop + 25]
        # once the exon is fully created we delete the gene sequence for memory efficiency
        printd("3' ss sequence")
        printd(self.seq_3ss)
        printd("5' ss sequence")
        printd(self.seq_5ss)
        self.gene.sequence = None


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

