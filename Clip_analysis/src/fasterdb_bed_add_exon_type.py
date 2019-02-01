#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    Add the exon type to the column name in a bed file containing the coordinates of fasterDB exons + add 200bp of \
    intron sequence surrounding the exons
"""

import sqlite3
import subprocess
import os


def get_exon(cnx):
    """
    Get every exons info to create a bed file
    :param cnx: (sqlite3 connect object) connection to fasterDB lite database
    :return: (list of string) list of the line tha will be display in the new bed
    """
    cursor = cnx.cursor()
    query = """SELECT t1.chromosome, t1.start_on_chromosome, t1.end_on_chromosome, t2.official_symbol, t1.pos_on_gene,
               t1.exon_type, t2.strand, t1.end_on_chromosome - t1.start_on_chromosome + 1
               FROM exons t1, genes t2
               WHERE t1.id_gene = t2.id
               ORDER BY t1.chromosome ASC, t1.start_on_chromosome ASC
            """
    cursor.execute(query)
    res = cursor.fetchall()
    list_res = []
    for exon in res:
        exon = list(exon)
        if int(exon[-1]) > 0:
            if not exon[5]:
                exon[5] = "NA"
            if "-" in str(exon[6]):
                exon[6] = "-"
            else:
                exon[6] = "+"
            new_exon = exon[0:3] + ["%s_%s_%s" % (exon[3], exon[4], exon[5])] + [".", exon[6], exon[7]]
            list_res.append("\t".join(list(map(str, new_exon))))
    return list_res


def get_intron_sequence(cnx):
    """
    Return every intron sequence in fasterDB
    :param cnx:  (sqlite3 connect object) connection to fasterDB database
    :return: (list of string) list of every intergenic sequence
    """
    cursor = cnx.cursor()
    query = """SELECT t2.chromosome, t1.start_on_chromosome,
                      t1.end_on_chromosome, t2.official_symbol, t1.pos_on_gene, t2.strand,
                      t1.end_on_gene - t1.start_on_gene + 1
               FROM introns t1, genes t2
               WHERE t1.id_gene = t2.id
               ORDER BY t2.chromosome ASC, t1.start_on_chromosome ASC"""
    cursor.execute(query)
    res = cursor.fetchall()
    print("introns : %s" % len(res))
    list_res = []
    for exon in res:
        exon = list(exon)
        if int(exon[-1]) > 0:
            if "-" in str(exon[5]):
                exon[5] = "-"
            else:
                exon[5] = "+"
            new_exon = exon[0:3] + ["%s_%s" % (exon[3], exon[4])] + [".", exon[5], exon[6]]
            list_res.append("\t".join(list(map(str, new_exon))))
    return list_res


def write_new_bed(cnx, bed_file_exon_type, type_file="exon"):
    """
    Create the new bed file
    :param cnx: (sqlite3 connect object) connection to fasterDB lite database
    :param bed_file_exon_type:  (string) the resulting bed file with exon type
    :param type_file: (string) the type of feature to write in the result file.
    """
    if type_file == "exon":
        list_ft = get_exon(cnx)
    else:
        list_ft = get_intron_sequence(cnx)
        print(len(list_ft))
    with open(bed_file_exon_type, "w") as bedout:
        for i in range(len(list_ft)):
            bedout.write("%s\n" % list_ft[i])


def main():
    """
    Create a bed file with exon type and 200bp of surrouding intron sequence and a bed files of intragenic sequence
    """
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/bed_template/")
    if not os.path.isdir(output):
        os.mkdir(output)
    fasterdb_file = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/fasterDB_lite.db")
    final_exon_bed_file = "%s/fasterDB_exons.bed" % output
    final_gene_bed_file = "%s/fasterDB_intron.bed" % output
    cnx = sqlite3.connect(fasterdb_file)
    print("Writing intron bed file")
    write_new_bed(cnx, final_exon_bed_file)
    print("Writing intron bed file")
    write_new_bed(cnx, final_gene_bed_file, "intron")


if __name__ == "__main__":
    main()
