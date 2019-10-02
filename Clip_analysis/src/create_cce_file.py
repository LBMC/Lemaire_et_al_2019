#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
Cretae a file containing every CCE exon not regulated by any splicing factors.
"""

import os
import sqlite3
import union_dataset_function as udf


def get_control_exon_information(cnx, exon_type, exon2remove):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of list of 2int) list of exon to remove frome the control list of exons
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   WHERE t1.exon_type LIKE '%{}%'
                   AND t1.id_gene = t2.id""".format(exon_type)
    else:
        query = """SELECT t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   AND t1.id_gene = t2.id
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    print("number of control exon before removing bad ones : %s" % len(result))
    nresult = ["_".join(list(map(str, exon))) for exon in result if [exon[0], exon[1]] not in exon2remove]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    return nresult



def main():
    base_dir = os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)))
    fasterdb = base_dir + "/data/fasterDB_lite.db"
    seddb = base_dir + "/data/sed.db"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    exon2remove = udf.get_exon_regulated_by_sf(cnx_sed, "down")
    my_exons = get_control_exon_information(cnx, "CCE", exon2remove)
    with open("data/input/CCE_exons.txt", "w") as outfile:
        outfile.write("\n".join(my_exons) + "\n")


if __name__ == "__main__":
    main()