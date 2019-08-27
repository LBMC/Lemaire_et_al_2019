#!/usr/bin/bash

"""
Description:
    Convert microexons supplementary file into an exon bed
"""


import sqlite3
import re
import os


def find_genes(cnx, chromosome, start, stop, name):
    """

    :param cnx: (sqlite3 connection object) connection to fasterdb database
    :param chromosome: (str) the chromosome
    :param start: (int) start position of the gene
    :param stop: (int)  stop position of the gene
    :param name: (str) the name of the gene
    """
    cursor = cnx.cursor()
    query = """
            SELECT strand, official_symbol
            FROM genes
            WHERE chromosome = '{0}'
            AND start_on_chromosome < {1}
            AND end_on_chromosome > {2}
    """.format(chromosome, start, stop)

    cursor.execute(query)
    res = cursor.fetchall()

    if len(res) == 0:
        print("Warning : %s chr%s:%s-%s doesn't have any gene related"
              % (name, chromosome, start, stop))
        return None
    elif len(res) == 1:
        mon_gene = res[0]
        if mon_gene[1] != name:
            print("Warning : %s chr%s:%s-%s have only one gene matching "
                  "with it's coordinates but with a different name : %s != %s"
                  % (name, chromosome, start, stop, name, mon_gene[1]))
            return None
        return mon_gene[0]
    else:
        for mon_gene in res:
            if mon_gene[1] == name:
                return mon_gene[0]
        print("Warning %s genes found matching %s chr%s:%s-%s "
              "but none have the good name" %
              (len(res), name, chromosome, start, stop))
        print(res)
        return res[0][0]


def write_bed(cnx, mfile, output):
    """
    Write a bed file of the micro-exons found by irimia et al.

    :param cnx: (sqlite3 connection object) connection to fasterdb database
    :param mfile: (str) the suplemental file of exons
    :param output: (str) folder where the result will be created
    """
    dics = {-1: "-", 1: "+"}
    dic_name = {}
    infile = open(mfile, "r")
    with open("%s/Irimia_et_al_microexons.bed" % output, "w") as outfile:
        for line in infile:
            if "#" not in line:
                line = line.replace("\n", "").split("\t")
                coordinates = re.split(r":|-", line[1].replace("chr", ""))
                chromosome = coordinates[0]
                start = int(coordinates[1]) - 1
                stop = int(coordinates[2])
                exon_name = line[0]
                i = 1
                while exon_name + "_" + str(i) in dic_name.keys():
                    i += 1
                dic_name[exon_name + "_" + str(i)] = 1
                strand = find_genes(cnx, chromosome, start, stop, line[0])
                if strand is not None:
                    outfile.write("\t".join([chromosome, str(start), str(stop),
                                             exon_name + "_" + str(i),
                                             "0", dics[strand]]) + "\n")

    infile.close()


def main():
    base = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    seddb = base + "/data/fasterDB_lite.db"
    output = base + "/result/irimia_bed/"
    mfile = base + "/data/irimia_et_al/table_S2.csv"
    if not os.path.isdir(output):
        os.mkdir(output)
    cnx = sqlite3.connect(seddb)
    write_bed(cnx, mfile, output)


if __name__ == "__main__":
    main()
