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


def get_intragenic_sequence(cnx):
    """
    Return every intragenic sequence in fasterDB
    :param cnx:  (sqlite3 connect object) connection to fasterDB database
    :return: (list of string) list of every intergenic sequence
    """
    cursor = cnx.cursor()
    query = """SELECT chromosome, start_on_chromosome, end_on_chromosome, official_symbol, strand, 
                      end_on_chromosome - start_on_chromosome + 1
               FROM genes
               ORDER BY chromosome ASC, start_on_chromosome ASC"""
    cursor.execute(query)
    res = cursor.fetchall()
    list_res = []
    for gene in res:
        gene = list(gene)
        if int(gene[-1]) > 0:
            if "-" in str(gene[4]):
                gene[4] = "-"
            else:
                gene[4] = "+"
            list_res.append("\t".join(list(map(str, gene))))
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
        list_ft = get_intragenic_sequence(cnx)
        print(len(list_ft))
    with open(bed_file_exon_type, "w") as bedout:
        for i in range(len(list_ft)):
            bedout.write("%s\n" % list_ft[i])


def add_intron_sequence(bed_file_exon_type, final_bed_file, chrom_size_file, left_intron=200, right_intron=200):
    """
    Add 200nt of intron sequence before and after each exons
    :param bed_file_exon_type: (string) bed file with exon type
    :param final_bed_file: (string) bed file with exon type + 200 nt of sourroungin intron sequence for each exon
    :param chrom_size_file: (string) a file containing chromosome size
    :param left_intron: (int) nt add to the left border of the intron
    :param right_intron: (int) nt add to the right location
    """
    cmd = "bedtools slop -i %s -g %s -l %s -r %s > %s" % (bed_file_exon_type, chrom_size_file, left_intron,
                                                          right_intron, final_bed_file)
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
    cmd = "rm %s" % bed_file_exon_type
    #subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)


def main():
    """
    Create a bed file with exon type and 200bp of surrouding intron sequence and a bed files of intragenic sequence
    """
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/bed_template/")
    if not os.path.isdir(output):
        os.mkdir(output)
    fasterdb_file = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/fasterDB_lite.db")
    chrom_size_file = fasterdb_file.replace("fasterDB_lite.db", "hg19.ren.chrom.sizes")
    final_exon_bed_file = "%s/fasterDB_exons_add200nt.bed" % output
    final_gene_bed_file = "%s/fasterDB_gene.bed" % output
    cnx = sqlite3.connect(fasterdb_file)
    bed_file_exon_type = "%s/tmp.bed" % output
    print("Writing intron bed file")
    write_new_bed(cnx, bed_file_exon_type)
    add_intron_sequence(bed_file_exon_type, final_exon_bed_file, chrom_size_file)
    print("Writing intragenic_bed file")
    write_new_bed(cnx, final_gene_bed_file, "gene")


if __name__ == "__main__":
    main()
