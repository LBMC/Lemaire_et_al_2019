#!/usr/bin/env python3


# -*- coding: utf-8 -*-

import os
import sqlite3
import exon_class_metaexon
import win_size
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)).replace("/metaexon_figure", ""))
import union_dataset_function


def get_control_exon_information(cnx, exon_type, exon2remove):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of 2 int) the exons we want to remove
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


def control_dictionaries_creator(window_size):
    """
    Create the control dictionary containing the vector of values of control exons. those vector will be use to \
    display the frequencies of a given nucleotide in a meta-exon figures for the control exons. \
    Create control dictionary files that contain the values for the boxplot, metagene and metagene windowed figure for \
    coding CCE exons, CE exon, ACE exons and ASE exons.
    :param window_size: (int) the size of the window we want to use to create the control metagene windowsed \
    dictionaries
    """
    exon_class_metaexon.set_debug(0)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fasterdb = os.path.dirname(os.path.realpath(__file__)).replace("src/metaexon_figure", "data/fasterDB_lite.db")
    seddb = os.path.dirname(os.path.realpath(__file__)).replace("src/metaexon_figure", "data/sed.db")
    ctrl_dir = dir_path + "/control_dictionaries/"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    if not os.path.isdir(ctrl_dir):
        os.mkdir(ctrl_dir)
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx_sed, "down")
    exon_type = ["ACE", "CCE"]
    for cur_exon_type in exon_type:
        ctrl_exon_list = get_control_exon_information(cnx, cur_exon_type, exon2remove)
        # ctrl_exon_list = ctrl_exon_list[0:2]
        list_exon = [exon_class_metaexon.ExonClass(cnx, exon[0], exon[1], exon[2], window_size) for exon in ctrl_exon_list]
        print("creating metagene windowsed")
        final_res_5p, final_res_3p, p5_analyzed, p3_analyzed = \
            exon_class_metaexon.get_metagene_vectors_windowsed(list_exon, window_size)
        print(final_res_5p)
        cur_file = open(ctrl_dir + cur_exon_type + "_metagene_windowsed.py", "w")
        cur_file.write("final_res_5p=" + str(final_res_5p) + "\n")
        cur_file.write("final_res_3p=" + str(final_res_3p) + "\n")
        cur_file.write("# " + str(p5_analyzed) + " sequences 5' analysees\n")
        cur_file.write("# " + str(p3_analyzed) + " sequences 3' analysees\n")
        cur_file.close()
        del (final_res_5p, final_res_3p, p5_analyzed, p3_analyzed)


def main():
    """
    launch the creation of control dictionaries
    """
    window_size = win_size.window_size
    control_dictionaries_creator(window_size)


if __name__ == "__main__":
    main()
