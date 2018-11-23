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


# def get_exon_type(cnx, coordinates):
#     """
#     From chromosomal coordinates give the type of the exons with the coordinates ``coordinates
#     :param cnx: (sqlite3 connect object) connection to fasterDB lite database
#     :param coordinates: (list of onj string and 2 int) chromosome number + chromosomal coordinates
#     :return: (string) the exon type of the exons with the ``cordinates``
#     """
#     cursor = cnx.cursor()
#     query = "SELECT exon_type FROM exons WHERE chromosome = '?' AND start_on_chromosome = ? AND end_on_chromosome = ?"
#     cursor.execute(query, (coordinates[0], coordinates[1], coordinates[2]))
#     res = cursor.fetchone()
#     if res[0]:
#         return res[0]
#     else:
#         return "NA"
#
#
# def write_new_bed(cnx, bed_file, bed_file_exon_type):
#     """
#     Create the new bed file
#     :param cnx: (sqlite3 connect object) connection to fasterDB lite database
#     :param bed_file: (string) a bed file containing fasterDB exons
#     :param bed_file_exon_type:  (string) the resultting bed file with exon type
#     """
#     with open(bed_file, "r") as bedin, open(bed_file_exon_type, "w") as bedout:
#         for line in bedin:
#             if "#" not in line:
#                 line = line.split("\t")
#                 coordinates = [line[0], int(line[1]), int(line[2])]
#                 exon_type = get_exon_type(cnx, coordinates)
#                 line[3] = line[3] + "_" + exon_type
#                 bedout.write("\t".join(line))
#             else:
#                 bedout.write(line)


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
                exon[5] =  "NA"
            if "-" in str(exon[6]):
                exon[6] = "-"
            else:
                exon[6] = "+"
            new_exon = exon[0:3] + ["%s_%s_%s" % (exon[3], exon[4], exon[5])] + [".", exon[6], exon[7]]
            list_res.append("\t".join(list(map(str, new_exon))))
    return list_res


def write_new_bed(cnx, bed_file_exon_type):
    """
    Create the new bed file
    :param cnx: (sqlite3 connect object) connection to fasterDB lite database
    :param bed_file_exon_type:  (string) the resulting bed file with exon type
    """
    list_exon = get_exon(cnx)
    with open(bed_file_exon_type, "w") as bedout:
        bedout.write("\n".join(list_exon))


def add_intron_sequence(bed_file_exon_type, final_bed_file, chrom_size_file):
    """
    Add 200nt of intron sequence before and after each exons
    :param bed_file_exon_type: (string) bed file with exon type
    :param final_bed_file: (string) bed file with exon type + 200 nt of sourroungin intron sequence for each exon
    :param chrom_size_file: (string) a file containing chromosome size
    """
    cmd = "bedtools slop -i %s -g %s -l 200 -r 200 > %s" % (bed_file_exon_type, chrom_size_file, final_bed_file)
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
    cmd = "rm %s" % bed_file_exon_type
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)


def main():
    """
    Create a bed file with exon type and 200bp of surrouding intron sequence
    """
    fasterdb_file = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/fasterDB_lite.db")
    chrom_size_file = fasterdb_file.replace("fasterDB_lite.db", "hg19.ren.chrom.sizes")
    final_bed_file = fasterdb_file.replace("fasterDB_lite.db", "fasterDB_exons_add200nt.bed")
    cnx = sqlite3.connect(fasterdb_file)
    bed_file_exon_type = fasterdb_file.replace("fasterDB_lite.db", "tmp.bed")
    write_new_bed(cnx, bed_file_exon_type)
    add_intron_sequence(bed_file_exon_type, final_bed_file, chrom_size_file)


if __name__ == "__main__":
    main()


