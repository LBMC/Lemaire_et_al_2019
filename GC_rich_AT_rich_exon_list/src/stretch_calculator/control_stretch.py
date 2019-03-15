#!/usr/bin/env python3


# -*- coding: utf-8 -*-

"""
Description:
    Create control files to store the ppt and bp score.
"""

import os
import sqlite3
import stretch_calculator
import sys
import config
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("stretch_calculator",
                                                                       "make_control_files_bp_ppt/"))
import exon_class
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)).replace("/stretch_calculator", ""))
import union_dataset_function


def get_control_exon_information(cnx, exon_type, exon2remove):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of list of 2 int) list of exon we want to remove
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
    print("Number of control exons before removing the one regulated by a splicing factors : %s" % len(result))
    nresult = [list(exon) for exon in result if [exon[1], exon[2]] not in exon2remove]
    print("Number of control exons after removing the one regulated by a splicing factors : %s" % len(nresult))
    return nresult


def control_dictionaries_creator():
    """
    Create the control dictionary containing the values corresponding to the score of bp and ppt for every control exons
    """
    exon_class.set_debug(0)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fasterdb = os.path.dirname(os.path.realpath(__file__)).replace("src/stretch_calculator", "data/fasterDB_lite.db")
    seddb = os.path.dirname(os.path.realpath(__file__)).replace("src/stretch_calculator", "data/sed.db")
    ctrl_dir = dir_path + "/control_dictionaries/"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx_sed, "down")
    if not os.path.isdir(ctrl_dir):
        os.mkdir(ctrl_dir)
    exon_type = config.exon_type
    stretches = config.stretches
    sequence_boundaries = config.sequence_boundaries
    for cur_exon_type in exon_type:
        print("Working on %s exons ..." % cur_exon_type)
        ctrl_exon_list = get_control_exon_information(cnx, cur_exon_type, exon2remove)
        print("\t--> Getting upstream intron sequence")
        list_exon = [exon_class.ExonClass(cnx, exon[0], exon[1], exon[2]) for exon in ctrl_exon_list]
        cur_file = open(ctrl_dir + cur_exon_type + "_stretches.py", "w")
        for cur_stretch in stretches:
            print("\t--> Calculating stretches %s/%s" % (cur_stretch[1], cur_stretch[0]))
            stretch_dic = stretch_calculator.stretch_counter(list_exon, cur_stretch, sequence_boundaries)
            cur_file.write("stretch_%sX%s = %s\n" % (cur_stretch[0], cur_stretch[1], str(stretch_dic)))
        cur_file.close()


def main():
    """
    launch the creation of control dictionaries
    """
    control_dictionaries_creator()


if __name__ == "__main__":
    main()
