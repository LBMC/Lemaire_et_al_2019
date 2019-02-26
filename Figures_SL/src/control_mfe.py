#!/usr/bin/env python3


# -*- coding: utf-8 -*-

"""
Description:
    Create control files to store the ppt and bp score.
"""

import os
import sqlite3
import exon_class_mfe
import control_bp_ppt
import union_dataset_function
import function_mfe


def control_dictionaries_creator():
    """
    Create the control dictionary containing the values corresponding to the score of bp and ppt for every control exons
    """
    exon_class_mfe.set_debug(0)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fasterdb = os.path.dirname(os.path.realpath(__file__)).replace("src", "data/fasterDB_lite.db")
    seddb = os.path.dirname(os.path.realpath(__file__)).replace("src", "data/sed.db")
    ctrl_dir = dir_path + "/control_dictionaries/"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    if not os.path.isdir(ctrl_dir):
        os.mkdir(ctrl_dir)
    exon_type = "CCE"
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx_sed, "down")
    ctrl_exon_list = control_bp_ppt.get_control_exon_information(cnx, exon_type, exon2remove, "down")
    print("retrieving upstream intron sequence")
    list_exon = [exon_class_mfe.ExonClass(cnx, exon[0], exon[1], exon[2]) for exon in ctrl_exon_list]
    print("calculating mfe")
    mfe_list_3ss,  mfe_list_5ss= function_mfe.mfe_calculator(list_exon)
    cur_file = open(ctrl_dir + exon_type + "_mfe.py", "w")
    cur_file.write("mfe_3ss=" + str(mfe_list_3ss) + "\n")
    cur_file.write("mfe_5ss=" + str(mfe_list_5ss) + "\n")
    cur_file.close()


def main():
    """
    launch the creation of control dictionaries
    """
    control_dictionaries_creator()


if __name__ == "__main__":
    main()
