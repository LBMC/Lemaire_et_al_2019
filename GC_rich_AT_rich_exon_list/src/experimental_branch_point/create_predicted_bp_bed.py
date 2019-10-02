#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
The goal of this script is to launch svm bp finder on every intronic \
sequences of AT-exons, GC-exons and CCE-exons.
"""


import os
import sys
import sqlite3
from GC_AT_analysis_experimental_bp import read_file, get_ctrl_exons
src_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, src_folder + "/make_control_files_bp_ppt")
import exon_class_bp
import function_bp
sys.path.insert(0, src_folder)
import union_dataset_function as udf


def main():
    regulation = "down"
    exon_class_bp.set_debug(0)
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    output = base + "/result/experimental_branch_point"
    at_exon_file = base + "/result/AT_rich_exons"
    gc_exon_file = base + "/result/GC_rich_exons"
    fasterdb = base + "/data/fasterDB_lite.db"
    seddb = base + "/data/sed.db"
    cnx = sqlite3.connect(fasterdb)
    cnx_sed = sqlite3.connect(seddb)
    exon_type = "CCE"
    at_exon = read_file(at_exon_file)
    gc_exon = read_file(gc_exon_file)
    exon2remove = [list(map(int, exon))
                   for exon in udf.get_exon_regulated_by_sf(cnx_sed,
                                                            regulation)]
    ctrl_exons = get_ctrl_exons(cnx, exon_type, exon2remove)
    exon_list = gc_exon + at_exon + ctrl_exons
    type_exon = ["GC-exons"] * len(gc_exon) + \
                ["AT-exons"] * len(at_exon) + \
                ["%s-exons" % exon_type] * len(ctrl_exons)
    tot = len(exon_list)
    count = 0
    count_none = 0
    print("Creating bed of predicted branch points")
    with open("%s/predicted_branch_points.bed" % output, "w") as outf:
        for exon, name_exon in zip(exon_list, type_exon):
            exon = exon_class_bp.ExonClass(cnx, str(exon[0]), exon[0], exon[1])
            nb_good_bp, list_pos = function_bp.goob_bp_only(exon)
            if list_pos is not None:
                for line in list_pos:
                    line[3] += "_" + name_exon
                    line[0] = "chr" + str(line[0])
                    outf.write("\t".join(list(map(str, line))) + "\n")
            else:
                count_none += 1
            count += 1
            sys.stdout.write("%s/%s  (%s)              \r" %
                             (count, tot, count_none))
    cnx.close()
    cnx_sed.close()


if __name__ == "__main__":
    main()
