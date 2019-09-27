#!/usr/bin/env python3.5

"""
Description:
    Create a bed file containing every GC and AT exons and their GC frequency
"""

import sqlite3
import os
import sys
mydir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, mydir)
import union_dataset_function as udf
from figure_creator import get_exons_list
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
    nresult = [list(exon) for exon in result if [exon[0], exon[1]]
               not in exon2remove]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    return nresult


def is_in(exon, list_exons):
    """
    Return 1 if ``exon`` is in ``list_exons`` 0 else.

    :param exons: (list of 2 int)
    :param list_exons: (list of list of 2 int) list of exons
    :return: (int) Return 1 if ``exon`` is in ``list_exons`` 0 else.
    """
    return 1 if exon in list_exons else 0


def get_exon_info(cnx, sedb, fasterdb_file, exon_list, u1_exons, u2_exons):
    """

    :param cnx: (sqlite3 connect object) connexion to fasterdb
    :param fasterdb_file: (str) an sqlite3 database file
    :param sedb: (str) path to sed database
    :param exon_list:  (list of 2 int) list of exons
    :param u1_exons: (list of list of 2 int) list of exons regulated by U1
    :param u2_exons: (list of list of 2 int) list of exons regulated by U2
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
            stretch = catch_index_error(stretch_counter([cexon])["T"], 0)
            dic_info = {"GC_content": exon[6].split(";")[4],
                        "upstream_intron_size": exon[7],
                        "downstream_intron_size": exon[8],
                        "UNA_count": catch_index_error(exon_data[8], 0),
                        "Hbound_count": catch_index_error(exon_data[6], 0),
                        "good_bp": catch_index_error(exon_data[3], 0),
                        "MFE_5SS": catch_index_error(mfe_5ss, 0),
                        "MFE_3SS": catch_index_error(mfe_3ss, 0),
                        "T_stretch": stretch,
                        "U1-regulated": is_in(exon[0:2], u1_exons),
                        "U2-regulated": is_in(exon[0:2], u2_exons),
                        }

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
                    "T_stretch": stretch,
                    "U1-regulated": is_in(exon[0:2], u1_exons),
                    "U2-regulated": is_in(exon[0:2], u2_exons),
                    }
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
    if not os.path.isdir(output):
        os.mkdir(output)
    cnx = sqlite3.connect(seddb)
    u1_exons = [list(map(int, exon)) for exon in
        get_exons_list(cnx, ["SNRPC", "SNRNP70", "DDX5_DDX17"], "down")]
    u2_exons = [list(map(int, exon)) for exon in
                get_exons_list(cnx, ["U2AF2", "SF1", "SF3A3", "SF3B4"], "down")]
    print("U1-exons : %s exons" % len(u1_exons))
    print("U2-exons : %s exons" % len(u2_exons))
    exon_list = udf.get_exon_regulated_by_sf(cnx, "down")
    print("Getting exon data ...")
    exon_data = get_exon_info(cnx, seddb, fasterdb, exon_list, u1_exons,
                              u2_exons)
    print("Writing bed")
    write_bed(output, exon_data, "data_for_regulated_exons")

    cnx.close()
    cnx = sqlite3.connect(seddb)
    exon2remove = udf.get_exon_regulated_by_sf(cnx, "down")
    cnx_fasterdb = sqlite3.connect(fasterdb)
    exon_list = get_control_exon_information(cnx_fasterdb, "CCE", exon2remove) + exon2remove
    cnx_fasterdb.close()
    print("CCE exons + regulated exons : %s" % len(exon_list))
    print("Getting CCE + regulated exon data ...")
    exon_data = get_exon_info(cnx, seddb, fasterdb, exon_list, u1_exons, u2_exons)
    print("Writing bed")
    write_bed(output, exon_data, "data_for_regulated_CCE_exons")
    cnx.close()


if __name__ == "__main__":
    main()