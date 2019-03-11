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
    def __init__(self, gene_name, gene_id, exon_position):
        """
        Initiate the creation of an exon

        :param gene_name: (string) official symbol for an exons
        :param gene_id: (int) the id of a gene
        :param exon_position:  (int) the position of the exon on the gene
        """
        self.gene = Gene(gene_name, gene_id)
        self.position = exon_position


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
        ExonClassMain.__init__(self, gene_name, gene_id, exon_position)
        self.gene.gene_filler(cnx)
        self.upstream_intron = Intron(cnx, gene_id, self.position - 1, "upstream", self.gene.sequence)
        # once the exon is fully created we delete the gene sequence for memory efficiency
        printd("Upstream proxi sequence")
        printd(self.upstream_intron.sequence_proxi)
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


class Intron:
    """Create an intron"""
    def __init__(self, cnx, gene_id, pos_on_gene, location, gene_seq):
        """

        :param cnx: (sqlite3 object) allow connection to **FasterDB Lite** database
        :param gene_id:  (int) the id of a gene
        :param pos_on_gene: (int) an intron position on a gene
        :param location: (string) upstream or downstream
        :param gene_seq: (string) the sequence of the gene
        """
        self.gene_id = gene_id
        self.position = pos_on_gene
        self.location = location
        self.start = None
        self.end = None
        self.sequence_proxi = None
        self.get_intron_information(cnx, gene_seq)

    def get_intron_information(self, cnx, gene_seq):
        """
        Get the intron coordinates on it's gene

        :param cnx: (sqlite3 object) allows connection to **FasterDB Lite** database
        :param gene_seq: (string) the sequence of the gene
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
            return None
        start = result[0][0] - 1
        end = result[0][1]
        if start < end:
            self.start = start
            self.end = end
        else:
            return None
        sequence = gene_seq[start:end]
        if len(sequence) > 19:
            if len(sequence) < 100:
                self.sequence_proxi = sequence
            else:
                self.sequence_proxi = sequence[-100:]


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
