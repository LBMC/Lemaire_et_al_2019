#!/usr/bin/env python3

"""
Description:
    The goal of this script is to compute the GC content of every exons \
    in a bed file (with hg19 coordinates) and get the length of \
    their flanking introns using fasterdb annotation.
"""

import sqlite3
from pyfaidx import Fasta
from Bio.Seq import Seq
import os
import lazyparser as lp


class NoIntersectionError(Exception):
    pass


class GeneNameError(Exception):
    pass


class AmbigousIntervalError(Exception):
    pass


def handle_query_result(res, namee, verbose, my_func="exon"):
    """
    Check the intersection result research.

    :param res: (list of list) the exons/introns intersecting the interval of \
    interest.
    :param namee: (str) the name of the exon
    :param verbose: (boolean) verbose mode
    :param my_func: (str) the name of the function where the result came from.
    :return: (list of list) the list of exons/intron intersection the \
    interval.
    """
    if len(res) == 0:
        if my_func == "exon":
            if verbose:
                print("\t\tNo intersection found")
            return None
        else:
            raise NoIntersectionError("No intersection found")
    elif len(res) == 1:
        interval = res[0]
        if namee.split("_")[0] != interval[1]:
            raise GeneNameError("Intersection Gene name is %s and the "
                                "bed gene name is %s" % (interval[1], namee))
        if verbose:
            print("\t\tone intersection found : %s" % interval)
        return interval
    else:
        if verbose:
            print("\t\t%s intersections found : %s" % (len(res), res))
        for interval in res:
            if namee.split("_")[0] == interval[1]:
                if verbose:
                    print("\t\tInterval %s selected" % interval)
                return interval
        raise AmbigousIntervalError("None of the intervals  %s have the gene "
                                    "name %s" % (res, namee))


def fasterdb_exon_intersection(cursor, chromosome, start, stop, namee, strand,
                               verbose=False):
    """
    Find the exons intersecting the interval defined by ``chromosome``,
    ``start`` and ``stop``.

    :param cursor: (sqlite3 cursor object) connection to fasterdb
    dic_strand
    :param chromosome: (str) the name of the chromosome of the interval \
    of interest (the interval composed of ``chromosome``, ``start``,
    ``stop`` should be coming from a bed file of exons.
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param namee: (str) the name of the exonx
    :param strand: (str) the strand
    :param verbose: (boolean) verbose mode
    :return: (list of list) the list of exons/intron intersection the \
    interval.
    """
    dic_strand = {"-": -1, "+": 1}
    query = """
            SELECT t1.id_gene, t2.official_symbol, t1.pos_on_gene, 
                   t2.chromosome, t1.start_on_chromosome - 1 , 
                   t1.end_on_chromosome,
                   t2.strand, t3.upstream_intron_size, 
                   t3.downstream_intron_size
            FROM exons t1, genes t2, s.sed t3
            WHERE t1.id_gene = t2.id
            AND t2.chromosome = '{0}'
            AND t2.strand = {1}
            AND t1.id_gene = t3.gene_id
            AND t1.pos_on_gene = t3.exon_pos
            AND ((t1.start_on_chromosome - 1 <= {2}
                 AND t1.end_on_chromosome > {2})
            OR (t1.start_on_chromosome - 1 < {3}
                 AND t1.end_on_chromosome >= {3})
            OR ({2} < t1.start_on_chromosome - 1 
                AND {3} > t1.end_on_chromosome))
    """.format(chromosome, dic_strand[strand], start, stop)

    if verbose:
        print("\t\tQuery to execute: \n\t%s" % query)

    cursor.execute(query)
    res = cursor.fetchall()
    res = [list(interval) for interval in res]
    return handle_query_result(res, namee, verbose)


def redefine_intron_size4exon(chromosome, start, stop, strand, interval,
                              verbose=False):
    """
    Adjust flanking intron size for an exon defined by \
    ``chromosome`` ``start`` ``stop`` and ``strand`` and an intersection \
    fasterdb exon defined by ``interval``.

    :param chromosome: (str) a chromosome
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param strand: (str) the strand
    :param interval: (list of value) interval intersection the exons \
    defined by the other parameters
    :param verbose: (str) verbose mode
    :return: (dictionary of 2 int) the size of the upstream flanking intron \
    and the size of the downstream flanking intron.
    """
    dic_strand = {"-": -1, "+": 1}
    if chromosome != interval[3]:
        raise NameError("chromosome %s and %s should be the same" %
                        (chromosome, interval[3]))
    if dic_strand[strand] != interval[6]:
        raise NameError("The strands %s and %s sould be the same" %
                        (dic_strand[strand], interval[6]))
    dic = {"upstream_intron_size": None, "downstream_intron_size": None}
    if interval[8] is not None:
        dic["downstream_intron_size"] = interval[8] + interval[5] - stop
        if dic["downstream_intron_size"] < 0:
            raise ValueError("downstream_intron_size should not have a "
                             "negative length (%s)" %
                             dic["downstream_intron_size"])
    if interval[7] is not None:
        dic["upstream_intron_size"] = interval[7] + start - interval[4]
        if dic["upstream_intron_size"] < 0:
            raise ValueError("upstream_intron_size should not have a "
                             "negative length (%s)" %
                             dic["upstream_intron_size"])
    if verbose:
        print("\t\tRedefining flanking intron size for exons %s:%s-%s "
              "intersecting with %s exons" % (chromosome, start, stop,
                                              interval))
        print("\t\tResult : %s" % dic)
    return dic


def fasterdb_intron_intersection(cursor, chromosome, start, stop, namee,
                                 strand, verbose=False):
    """
    Find the intron intersecting the interval defined by ``chromosome``,
    ``start`` and ``stop``.

    :param cursor: (sqlite3 cursor object) connection to fasterdb
    dic_strand
    :param chromosome: (str) the name of the chromosome of the interval \
    of interest (the interval composed of ``chromosome``, ``start``,
    ``stop`` should be coming from a bed file of exons.
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param namee: (str) the name of the exonx
    :param strand: (str) the strand
    :param verbose: (boolean) verbose mode
    :return: (list of tuple) the list of exons/intron intersection the \
    interval.
    """
    dic_strand = {"-": -1, "+": 1}
    query = """
            SELECT t1.id_gene, t2.official_symbol, t1.pos_on_gene, 
                   t2.chromosome, t1.start_on_chromosome - 1 , 
                   t1.end_on_chromosome, t2.strand
            FROM f.introns t1, f.genes t2
            WHERE t1.id_gene = t2.id
            AND t2.chromosome = '{0}'
            AND t2.strand = {1}
            AND ((t1.start_on_chromosome - 1 <= {2}
                 AND t1.end_on_chromosome > {2})
            AND (t1.start_on_chromosome - 1 < {3}
                 AND t1.end_on_chromosome >= {3}))
    """.format(chromosome, dic_strand[strand], start, stop)

    if verbose:
        print("\t\tQuery to execute: \n\t%s" % query)

    cursor.execute(query)
    res = cursor.fetchall()
    res = [list(interval) for interval in res]
    return handle_query_result(res, namee, verbose, my_func="intron")


def redefine_intron_size4intron(chromosome, start, stop, strand, interval,
                                verbose=False):
    """
    Adjust flanking intron size for an exon defined by \
    ``chromosome`` ``start`` ``stop`` and ``strand`` and an intersection \
    fasterdb intron defined by ``interval``.

    :param chromosome: (str) a chromosome
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param strand: (str) the strand
    :param interval: (list of value) interval intersection the exons \
    defined by the other parameters
    :param verbose: (boolean) verbose mode
    :return: (dictionary of 2 int) the size of the upstream flanking intron \
    and the size of the downstream flanking intron.
    """
    dic_strand = {"-": -1, "+": 1}
    if chromosome != interval[3]:
        raise NameError("chromosome %s and %s should be the same" %
                        (chromosome, interval[3]))
    if dic_strand[strand] != interval[6]:
        raise NameError("The strands %s and %s sould be the same" %
                        (dic_strand[strand], interval[6]))
    dic = {"upstream_intron_size": None, "downstream_intron_size": None}
    if interval[5] - stop > 0:
        dic["downstream_intron_size"] = interval[5] - stop
    elif interval[5] - stop < 0:
        raise ValueError("interval[5](%s) - stop(%s) should not be below "
                         "Zero" % (interval[5], stop))
    if start - interval[4] > 0:
        dic["upstream_intron_size"] = start - interval[4]
    elif start - interval[4] < 0:
        raise ValueError("start(%s) - interval[4](%s) should not be below "
                         "Zero" % (start, interval[4]))
    if verbose:
        print("\t\tRedefining flanking intron size for exons %s:%s-%s "
              "intersecting with %s introns" % (chromosome, start, stop,
                                                interval))
        print("\t\tResult : %s" % dic)
    return dic


def get_gc_content(chromosome, start, stop, strand, dic_seq):
    """
    Compute the GC content of the interval given by the first 4 parameters.

    :param chromosome: (str) a chromosome
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param strand: (str) the strand
    :param dic_seq: (pyfaidx dictionary) the human genome
    :return: (float) the gc content of the exons
    """
    sequence = Seq(str(dic_seq[chromosome][start:stop]))
    if strand == "-":
        sequence = sequence.reverse_complement()
    return round((sequence.count("G") + sequence.count("C"))
                 / len(sequence) * 100, 1)


def compute_exon_data(chromosome, start, stop, strand, namee, cursor, dic_seq,
                      verbose):
    """
    compute the upstream and downstream intron size and the gc content of \
    the exon.

    :param chromosome: (str) a chromosome
    :param start: (int) the start point of the interval on the ``chromosome``
    :param stop: (int) the stop point of the interval on ``chromosome``
    :param strand: (str) the strand
    :param namee: (str) the name of the exons
    :param dic_seq: (pyfaidx dictionary) the human genome
    :param verbose: (boolean) verbose mode
    :param cursor: (sqlite3 cursor) connection to fasterdb
    :return: (dictionary of int/float) dictionary that contains the \
    upstream and downstream intron size and the gc content of \
    the exon defined by the 5 first parameters
    """
    exon_name = "%s %s:%s-%s:%s" % (namee, chromosome, start, stop, strand)
    if verbose:
        print("Working on %s" % exon_name)
        print("\tChecking if exon %s overlaps an exon" % exon_name)
    interval = fasterdb_exon_intersection(cursor, chromosome, start, stop,
                                          namee, strand, verbose)
    if interval is None:
        if verbose:
            print("\tCheking if exons %s overlaps an intron" % exon_name)
        interval = fasterdb_intron_intersection(cursor, chromosome, start,
                                                stop,  namee, strand, verbose)
        if verbose:
            print("\tComputing flanking introns size")
        dic = redefine_intron_size4intron(chromosome, start, stop, strand,
                                          interval, verbose)
    else:
        if verbose:
            print("\tComputing flanking introns size")
        dic = redefine_intron_size4exon(chromosome, start, stop, strand,
                                        interval, verbose)
    if verbose:
        print("\tComputing GC content of exon %s" % exon_name)
    gc_content = get_gc_content(chromosome, start, stop, strand, dic_seq)
    dic["GC_content"] = gc_content
    if verbose:
        print("\tFinal data recovered for exons %s : %s" % (exon_name, dic))
    return dic


def write_new_bed(input_bed, output_file, cursor, dic_seq, verbose):
    """
    Write a new bed file with some exon data in a new 7th columns.

    :param input_bed: (str) the input be file of exons for which we want \
    to get exons data.
    :param output_file: (str) new bed file some exon data in a new 7th columns
    :param cursor: (sqlite3 cursor object) connection to fasterdb
    :param dic_seq: (pyfaidx dictionary) the human genome
    :param verbose: (boolean) verbose mode
    :return:
    """
    infile = open(input_bed, "r")
    with open(output_file, "w") as ouf:
        for line in infile:
            if "#" not in line:
                line = line.replace("\n", "")
                line = line.split("\t")
                data = compute_exon_data(line[0], int(line[1]), int(line[2]),
                                         line[5], line[3], cursor, dic_seq,
                                         verbose)
                line.append(str(data))
                ouf.write("\t".join(list(map(str, line))) + "\n")
    infile.close()


@lp.parse(input_bed="file", output="dir", faster_db="file", seddb="file",
          human_genome="file", verbose=["y", "n"])
def main(input_bed, output, faster_db, seddb, human_genome, verbose="n"):
    """
    Compute the GC content of every exons \
    in a bed file (with hg19 coordinates) and get the length of \
    their flanking introns using fasterdb annotation.

    :param input_bed: (str) The input bed containing micro-exons
    :param output: (str) Folder where the result will be created
    :param faster_db: (str) Path to fasterdb database
    :param seddb: (str) Path to sed database
    :param human_genome: (str) Human genome assembly fasta
    :param verbose: (str) y to enable verbose mode
    """
    vdic = {"y": True, "n": False}
    verbose = vdic[verbose]
    cnx = sqlite3.connect(faster_db)
    cursor = cnx.cursor()
    cursor.execute("ATTACH DATABASE ? as s", (seddb,))
    cursor.execute("ATTACH DATABASE ? as f", (faster_db,))
    dic_seq = Fasta(human_genome)
    output_file = "%s/%s_freq.bed" % (output, os.path.basename(input_bed)
                                      .replace(".bed", ""))
    write_new_bed(input_bed, output_file, cursor, dic_seq, verbose)
    cursor.close()
    cnx.close()


if __name__ == "__main__":
    main()
