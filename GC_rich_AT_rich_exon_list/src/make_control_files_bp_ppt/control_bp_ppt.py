#!/usr/bin/env python3


# -*- coding: utf-8 -*-

"""
Description:
    Create control files to store the ppt and bp score.
"""

import os
import sqlite3
import exon_class_bp
import function
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)).replace("/make_control_files_bp_ppt", ""))
import union_dataset_function


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
        query = """SELECT t2.official_symbol, t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   WHERE t1.exon_type LIKE '%{}%'
                   AND t1.id_gene = t2.id""".format(exon_type)
    else:
        query = """SELECT t2.official_symbol, t1.id_gene, t1.pos_on_gene
                   FROM exons t1, genes t2
                   AND t1.id_gene = t2.id
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    print("number of control exon before removing bad ones : %s" % len(result))
    nresult = [list(exon) for exon in result if [exon[1], exon[2]] not in exon2remove]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    return nresult


def control_dictionaries_creator():
    """
    Create the control dictionary containing the values corresponding to the score of bp and ppt for every control exons
    """
    exon_class.set_debug(0)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fasterdb = os.path.dirname(os.path.realpath(__file__)).replace("src/make_control_files_bp_ppt",
                                                                   "data/fasterDB_lite.db")
    seddb = os.path.dirname(os.path.realpath(__file__)).replace("src/make_control_files_bp_ppt", "data/sed.db")
    ctrl_dir = dir_path + "/control_dictionaries/"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx_sed, "down")
    if not os.path.isdir(ctrl_dir):
        os.mkdir(ctrl_dir)
    exon_type = ["CCE"]
    sizes = [100, 50, 35, 25]
    for cur_exon_type in exon_type:
        ctrl_exon_list = get_control_exon_information(cnx, cur_exon_type, exon2remove)
        print("retrieving upstream intron sequence")
        list_exon = [exon_class.ExonClass(cnx, exon[0], exon[1], exon[2]) for exon in ctrl_exon_list]
        for size in sizes:
            print("calculating bp and ppt score")
            bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, sequence_list, ag_count_list, \
                hbound_list, uaa_list, una_list = function.bp_ppt_calculator(list_exon, size)
            cur_file = open(ctrl_dir + cur_exon_type + "_" + str(size) + "_bp_ppt_score.py", "w")
            cur_file.write("bp_score=" + str(bp_score_list) + "\n")
            cur_file.write("ppt_score=" + str(ppt_score_list) + "\n")
            cur_file.write("nb_bp=" + str(nb_bp_list) + "\n")
            cur_file.write("nb_good_bp=" + str(nb_good_bp_list) + "\n")
            cur_file.write("bp_seq=" + str(sequence_list) + "\n")
            cur_file.write("ag_count=" + str(ag_count_list) + "\n")
            cur_file.write("hbound=" + str(hbound_list) + "\n")
            cur_file.write("uaa_count=" + str(uaa_list) + "\n")
            cur_file.write("una_count=" + str(una_list) + "\n")
            cur_file.close()


def main():
    """
    launch the creation of control dictionaries
    """
    control_dictionaries_creator()


if __name__ == "__main__":
    main()
