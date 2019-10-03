#!/usr/bin/env python3

"""
Description:
    The goal of this script is to create from the SED database to list of \
    exons : one corresponding to small exons and one corresponding to control \
    exons.
"""

import sqlite3
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import union_dataset_function as udf


dics = {-1: "-", 1: "+"}

def get_small_exons(cnx, seddb, fasterdb, threshold, exon_type, exon2remove):
    """
    Return every small exons data.

    :param cnx: (sqlite3 connection object) connection to sed database.
    :param seddb: (str) path to the simple exons description database
    :param fasterdb: (str) path to fasterdb
    :param threshold: (int) the size threshold
    :param exon_type: (str) the type of exon used
    :param exon2remove: (list of 2int) list of exon2remove
    :return: (list of small exons)
    """
    cursor = cnx.cursor()
    cursor.execute("ATTACH DATABASE ? as sed", (seddb,))
    cursor.execute("ATTACH DATABASE ? as fasterdb", (fasterdb,))

    query = """
           SELECT t1.id_gene, t1.pos_on_gene, t1.chromosome, 
                  t1.start_on_chromosome - 1, t1.end_on_chromosome,
                  t2.strand, t3.iupac_exon, t3.upstream_intron_size, 
                  t3.downstream_intron_size
           FROM fasterdb.exons t1, fasterdb.genes t2, sed.sed t3
           WHERE t1.id_gene = t2.id
           AND t1.id_gene = t3.gene_id
           AND t1.pos_on_gene = t3.exon_pos
           AND t3.exon_size <= %s 
           AND t3.exon_size > 2
    """ % threshold

    if exon_type == "ALL":
        query += ";"
    else:
        query += "AND t3.exon_type LIKE '%{0}%';".format(exon_type)

    cursor.execute(query)
    result1 = cursor.fetchall()
    result = [list(exon) for exon in result1]

    list_small_exons = [exon[2:5] + ["%s_%s" % (exon[0], exon[1]), ".",
                                     dics[exon[5]]]
                        + [str({"GC_content": float(exon[6].split(";")[4]),
                                "upstream_intron_size": exon[7],
                                "downstream_intron_size": exon[8]})]
                        for exon in result
                        if [exon[0], exon[1]] not in exon2remove]

    return list_small_exons


def get_small_sf_down_exons(cnx, seddb, fasterdb, threshold):
    """
    Return every small exons data.

    :param cnx: (sqlite3 connection object) connection to sed database.
    :param seddb: (str) path to the simple exons description database
    :param fasterdb: (str) path to fasterdb
    :param threshold: (int) the size threshold
    :return: (list of small exons)
    """
    cursor = cnx.cursor()
    # cursor.execute("ATTACH DATABASE ? as sed", (seddb,))
    # cursor.execute("ATTACH DATABASE ? as fasterdb", (fasterdb,))

    query = """
           SELECT t1.id_gene, t1.pos_on_gene, t1.chromosome, 
                  t1.start_on_chromosome - 1, t1.end_on_chromosome,
                  t2.strand, t3.iupac_exon, t3.upstream_intron_size, 
                  t3.downstream_intron_size
           FROM fasterdb.exons t1, fasterdb.genes t2, sed.sed t3
           WHERE t1.id_gene = t2.id
           AND t1.id_gene = t3.gene_id
           AND t1.pos_on_gene = t3.exon_pos
           AND t3.exon_size <= %s 
           AND t3.exon_size > 2;
    """ % threshold

    cursor.execute(query)
    result1 = cursor.fetchall()
    result = [list(exon) for exon in result1]

    down_sf_exons = udf.get_exon_regulated(cnx, "down")
    list_small_exons = [exon[2:5] + ["%s_%s" % (exon[0], exon[1]), ".",
                                     dics[exon[5]]]
                        + [str({"GC_content": float(exon[6].split(";")[4]),
                                "upstream_intron_size": exon[7],
                                "downstream_intron_size": exon[8]})] for exon
                        in result if [exon[0], exon[1]] in down_sf_exons]

    return list_small_exons


def get_control_exons(cnx, seddb, fasterdb, threshold, exon2remove, exon_type):
    """
    Return every control exons.

    :param cnx: (sqlite3 connection object) connection to sed database.
    :param seddb: (str) path to the simple exons description database
    :param fasterdb: (str) path to fasterdb
    :param threshold: (int) the size threshold
    :param exon2remove: (list of 2int) list of exon2remove
    :param exon_type: (str) the control exon type
    :return: (list of small exons)
    """

    cursor = cnx.cursor()
    # cursor.execute("ATTACH DATABASE ? as sed", (seddb,))
    # cursor.execute("ATTACH DATABASE ? as fasterdb", (fasterdb,))

    query = """
           SELECT t1.id_gene, t1.pos_on_gene, t1.chromosome, 
                  t1.start_on_chromosome - 1, t1.end_on_chromosome,
                  t2.strand, t3.iupac_exon, t3.upstream_intron_size, 
                  t3.downstream_intron_size
           FROM fasterdb.exons t1, fasterdb.genes t2, sed.sed t3
           WHERE t1.id_gene = t2.id
           AND t1.id_gene = t3.gene_id
           AND t1.pos_on_gene = t3.exon_pos
           AND t3.exon_size > {0}
        """.format(threshold)

    if exon_type == "ALL":
        query += ";"
    else:
        query += "AND t3.exon_type LIKE '%{0}%';".format(exon_type)

    cursor.execute(query)
    result1 = cursor.fetchall()
    result = [list(exon) for exon in result1]

    list_exons = [exon[2:5] + ["%s_%s" % (exon[0], exon[1]), ".", dics[exon[5]]]
                  + [str({"GC_content" : float(exon[6].split(";")[4]),
                      "upstream_intron_size": exon[7],
                       "downstream_intron_size": exon[8]})] for exon in
                  result if [exon[0], exon[1]] not in exon2remove]

    return list_exons


def write_bed_file(exon_list, output, name_fig):
    """
    Write a bed file.

    :param exon_list: (list of list of values) list of exons
    :param output: (str) path where the result will be created
    :param name_fig: (str) the name of the figure
    """
    with open("%s/%s" % (output, name_fig), "w") as outfile:
        for exon in exon_list:
            exon = list(map(str, exon))
            outfile.write("\t".join(exon) + "\n")


def main():
    """
    Create the graphics wanted
    """
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    seddb = base + "/data/sed.db"
    fasterdb = base + "/data/fasterDB_lite.db"
    output = base + "/result/variance_analysis/bed_file/"
    if not os.path.isdir(output):
        os.mkdir(output)
    exon_type = "CCE"
    regulation = "down"
    threshold = 27
    cnx = sqlite3.connect(seddb)
    exon2remove = udf.get_exon_regulated_by_sf(cnx, regulation)
    print("Get every FasterDB small exons : (exons with a size greater than 2"
          " nucleotides and lower than %s nucleotides" % threshold)
    small_exons = get_small_exons(cnx, seddb, fasterdb, threshold, exon_type,
                                  exon2remove)
    print("\tSmall exons found : %s" % len(small_exons))
    print("Getting small exons down-regulated by a splicing factor ...")
    r = [exon[3] for exon in small_exons]
    with open("%s/small.txt" % output, "w") as ouf:
        ouf.write("\n".join(r) + "\n")
    small_down_sf_exon = get_small_sf_down_exons(cnx, seddb, fasterdb,
                                                 threshold)
    print("\tSmall exons downregulated by a splicing factor found : %s"
          % len(small_down_sf_exon))
    print("Getting control exons ...")
    ctrl_exon = get_control_exons(cnx, seddb, fasterdb, threshold,
                                  exon2remove, exon_type)
    print("\tControl exons found : %s" % len(ctrl_exon))
    cnx.close()
    print("Writing results")
    write_bed_file(small_exons, output, "small_%s_exons_(3-%snt).bed" %
                   (exon_type, threshold))
    write_bed_file(small_down_sf_exon, output,
                   "small_sf-downregulated_exons_(3-%snt).bed" % threshold)
    write_bed_file(ctrl_exon, output, "%s_exons_%snt+.bed" %
                   (exon_type, threshold))


if __name__ == "__main__":
    main()