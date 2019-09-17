#!/usr/bin/env python3.5

"""
Description:
    Create a bed file containing every GC and AT exons and their GC frequency
"""

import sqlite3
import os
import numpy as np
import sys
import math
mydir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, mydir)
import union_dataset_function as udf
bp_dir = mydir + "/make_control_files_bp_ppt"
sys.path.insert(0, bp_dir)
import exon_class_bp
exon_class_bp.set_debug(0)
from function_bp import bp_ppt_calculator
mfe_dir = mydir + "/minimum_free_energy"
sys.path.insert(0, mfe_dir)
import exon_class
exon_class.set_debug(0)
from function import mfe_calculator
stretch_dir = mydir + "/stretch_calculator"
sys.path.insert(0, stretch_dir)
from stretch_calculator import stretch_counter




def get_exon_from_file(exon_file):
    """
    Get evey exon stored into ``exon_file``.

    :param exon_file: (str) an exon file
    :return: (list of 2 int) list of exons
    """
    exon_list = []
    with open(exon_file, "r") as infile:
        for line in infile:
            if line[0] != "#":
                line = line.replace("\n", "")
                line = line.split('\t')
                exon_list.append(list(map(int, line)))
    return exon_list


def catch_index_error(mlist, index):
    try:
        val = mlist[index]
    except IndexError:
        val = None
    return val

def get_exon_info(cnx, sedb, fasterdb_file, exon_list):
    """

    :param cnx: (sqlite3 connect object) connexion to fasterdb
    :param fasterdb_file: (str) an sqlite3 database file
    :param sedb: (str) path to sed database
    :param exon_list:  (list of 2 int) list of exons
    :return: (list of list of value) list of data
    """
    dic = {-1: "-", 1: "+"}
    cursor = cnx.cursor()
    cursor.execute("ATTACH DATABASE ? as sed", (sedb,))
    cursor.execute("ATTACH DATABASE ? as fasterdb", (fasterdb_file,))
    if exon_list is None:
        query = """
                SELECT t1.id_gene, t1.pos_on_gene, t1.chromosome, 
                       t1.start_on_chromosome, t1.end_on_chromosome, 
                       t2.strand, t3.iupac_exon, t3.upstream_intron_size,
                       t3.downstream_intron_size
                FROM fasterdb.exons as t1, fasterdb.genes as t2, sed.sed as t3
                WHERE t3.gene_id =  t1.id_gene
                AND   t3.exon_pos = t1.pos_on_gene
                AND   t1.id_gene = t2.id
                AND   t3.exon_type LIKE '%CCE%'
                """
        cursor.execute(query)
        res = cursor.fetchall()
        new_res = []
        for exon in res:
            exon = list(exon)
            exon[3] = int(exon[3]) - 1
            cexon = exon_class_bp.ExonClass(cnx, str(exon[0]), exon[0], exon[1])
            exon_data = bp_ppt_calculator([cexon])
            mexon = exon_class.ExonClass(cnx, str(exon[0]), exon[0], exon[1])
            mfe_5ss, mfe_3ss = mfe_calculator([mexon])
            stretch = catch_index_error(stretch_counter([cexon], [-75, -35])["T"], 0)
            dic_info = {"GC_content": exon[6].split(";")[4],
                        "upstream_intron_size": exon[7],
                        "downstream_intron_size": exon[8],
                        "UNA_count": catch_index_error(exon_data[8], 0),
                        "Hbound_count": catch_index_error(exon_data[6], 0),
                        "good_bp": catch_index_error(exon_data[3], 0),
                        "MFE_5SS": catch_index_error(mfe_5ss, 0),
                        "MFE_3SS": catch_index_error(mfe_3ss, 0),
                        "T_stretch": stretch}
            new_res.append(exon[2:5] + ["%s_%s" % (exon[0], exon[1])] + \
                    ["0", dic[exon[5]]] + [str(dic_info)])
        return new_res
    count = 0
    tot = len(exon_list)
    result = []
    for exon in exon_list:
        count += 1
        query = """
                SELECT t1.chromosome, t1.start_on_chromosome, 
                       t1.end_on_chromosome, t2.strand, t3.iupac_exon,
                       t3.upstream_intron_size, t3.downstream_intron_size
                FROM fasterdb.exons as t1, fasterdb.genes as t2, sed.sed as t3
                WHERE t3.gene_id =  t1.id_gene
                AND   t3.exon_pos = t1.pos_on_gene
                AND   t1.id_gene = t2.id
                AND   t3.gene_id = %s
                AND   t3.exon_pos = %s
                """ % (exon[0], exon[1])
        cursor.execute(query)
        res = cursor.fetchall()
        if len(res) > 1:
            raise IndexError("Error only one row shoud be return for %s" %
                             exon)
        tmp = list(res[0])
        tmp[1] = int(tmp[1]) - 1
        cexon = exon_class_bp.ExonClass(cnx, str(exon[0]), exon[0], exon[1])
        exon_data = bp_ppt_calculator([cexon])
        mexon = exon_class.ExonClass(cnx, str(exon[0]), exon[0], exon[1])
        mfe_5ss, mfe_3ss = mfe_calculator([mexon])
        stretch = catch_index_error(stretch_counter([cexon])["T"], 0)
        dic_info = {"GC_content": tmp[4].split(";")[4],
                    "upstream_intron_size": tmp[5],
                    "downstream_intron_size": tmp[6],
                    "UNA_count": catch_index_error(exon_data[8], 0),
                    "Hbound_count": catch_index_error(exon_data[6], 0),
                    "good_bp": catch_index_error(exon_data[3], 0),
                    "MFE_5SS": catch_index_error(mfe_5ss, 0),
                    "MFE_3SS": catch_index_error(mfe_3ss, 0),
                    "T_stretch": stretch}
        exon_data = tmp[0:3] + ["%s_%s" % (exon[0], exon[1])] + \
                    ["0", dic[tmp[3]]] + [str(dic_info)]
        result.append(exon_data)
        sys.stdout.write("Processing %s/%s\t\t\t\r" % (count, tot))
        sys.stdout.flush()
    return result


def write_bed(output, exon_data, name_bed):
    """
    Write a bed containing all the data present in ``exon_data``.

    :param output: (str) the folder were the result will be created
    :param exon_data: (list of list of value) list of data
    :param name_bed: (str) the name of the bed
    """
    with open("%s/%s.bed" % (output, name_bed), "w") as outfile:
        outfile.write("#chr\tstart\tstop\texon\tstrand\tscore\tGC_content\n")
        for exon in exon_data:
            outfile.write("\t".join(map(str, exon)) + "\n")


def main():
    """
    Create a bed file containing info about GC frequency of every GC-AT exons.
    """
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.realpath(__file__))))
    seddb = base + "/data/sed.db"
    fasterdb = base + "/data/fasterDB_lite.db"
    output = base + "/result/correlation_GC-AT-exons_TAD"
    at_rich_file = base + "/result/AT_rich_exons"
    gc_rich_file = base + "/result/GC_rich_exons"
    if not os.path.isdir(output):
        os.mkdir(output)
    cnx = sqlite3.connect(seddb)
    exon2remove = udf.get_exon_regulated_by_sf(cnx, "down")
    exons_at = get_exon_from_file(at_rich_file)
    print("Exons AT : %s" % len(exons_at))
    exons_gc = get_exon_from_file(gc_rich_file)
    print("Exons GC : %s" % len(exons_gc))
    exon_list = exons_at + exons_gc
    print("Getting exon data ...")
    exon_data = get_exon_info(cnx, seddb, fasterdb, exon_list)
    print("Writing bed")
    write_bed(output, exon_data, "GC_content_of_GC-AT_exons")

    cnx.close()
    cnx = sqlite3.connect(seddb)
    print("Getting CCE exon data ...")
    exon_data = get_exon_info(cnx, seddb, fasterdb, None)
    exon_list = ["%s_%s" % (exon[0], exon[1]) for exon in exon2remove]
    new_exon_data = [exon for exon in exon_data if exon[3] not in exon_list]
    print("Writing bed")
    write_bed(output, new_exon_data, "GC_content_of_CCE_exons")


if __name__ == "__main__":
    main()