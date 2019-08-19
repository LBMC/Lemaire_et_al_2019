#!/usr/bin env python3.5

"""
Description:
    The goal of this script is to create custom list of exons:
        - GA exons vs CT exons
        - GC-exons (small introns) not regulated by a splicing factor vs
          AT-exons (large introns) not regulated by a splicing
        - SRSF2-SRSF3-HNRNPC vs GC-exons (already used)
"""

import sqlite3
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import group_factor
from  figure_creator import get_exons_list


def file_writer(list_exon, name_list, output):
    """
    Write a file containing exons (on per line) identified by their gene_id \
    and their position within the hosting gene.

    :param list_exon: (list of 2 int) list of exons identified by their \
    gene_id and exon_position \
    within the hosting gene
    :param name_list: (string) the name of the list of exons named `list_exon`
    :param output: (string) path where the file will be created
    """
    with open("%s%s_exons.txt" % (output, name_list), "w") as my_file:
        for exon in list_exon:
            my_file.write("%s\t%s\n" % (exon[0], exon[1]))


def create_ct_ga_rich_exon_list(cnx, output):
    """
    Create the GA and the CT rich exons list.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param output: (str) path were the exon files will be created
    """
    ga_exon_all = get_exons_list(cnx, group_factor.ga_rich_down, "down")
    ct_exon_all = get_exons_list(cnx, group_factor.ct_rich_down, "down")
    ga_exon = [exon for exon in ga_exon_all if exon not in ct_exon_all]
    ct_exon = [exon for exon in ct_exon_all if exon not in ga_exon_all]
    print("ga_exon_all : %s exons" % len(ga_exon_all))
    print("ga_exon : %s exons" % len(ga_exon))
    print("ct_exon_all : %s exons" % len(ct_exon_all))
    print("ct_exon : %s exons" % len(ct_exon))
    file_writer(ga_exon, "GA_rich", output)
    file_writer(ct_exon, "CT_rich", output)


def venn_diagram_creator(list_exons, list_names, output):
    """
    Create a venn diagram of the `list1` and `list2` lists of values.

    :param list_exons: (list of list of 2 int) list of exons
    :param list_names: (list of string) the names of the list of exons
    :param output: (string) path where the result venn diagram will be created
    :return:
    """
    def flaten_list(x): return ["%s_%s" % (a[0], a[1]) for a in x]
    new_lists = [set(flaten_list(my_list)) for my_list in list_exons]
    venn3(new_lists, set_labels=list_names)
    plt.savefig("%sVenn_%s.pdf" % (output, "_".join(list_names)))
    plt.clf()
    plt.cla()
    plt.close()


def create_othergc_exon_file(cnx, output):
    """
    Create the GC-exons file regulated by SRSF2 HNRNPC and SRSF3.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param output: (str) path were the exon files will be created
    """
    gc_exon_all = get_exons_list(cnx, group_factor.gc_rich_down, "down")
    print("GC exons : %s" % len(gc_exon_all))
    other_gc_all = get_exons_list(cnx, group_factor.other, "down")
    print("Number of other_gc all exons : %s" % len(other_gc_all))
    other_gc = [exon for exon in other_gc_all if exon not in gc_exon_all]
    gc_exon = [exon for exon in gc_exon_all if exon not in other_gc_all]
    print("other GC exons : %s" % len(other_gc))
    file_writer(other_gc, "other_GC_rich", output)
    print("other GC exons all: %s" % len(other_gc_all))
    file_writer(other_gc_all, "other_GC_all_rich", output)
    print("GC exons: %s" % len(gc_exon))
    file_writer(gc_exon, "GC_rich", output)
    at_exon_all = get_exons_list(cnx, group_factor.at_rich_down, "down")
    print("Number of at exons all: %s" % len(at_exon_all))
    other_gc2 = [exon for exon in other_gc_all if exon not in at_exon_all]
    at_exon = [exon for exon in at_exon_all if exon not in other_gc_all]
    print("other GC exons (for AT comparison): %s" % len(other_gc2))
    print("Number of at exons : %s" % len(at_exon))
    file_writer(other_gc2, "other_GC_rich4ATcomp", output)
    file_writer(at_exon, "AT_rich", output)
    venn_diagram_creator([other_gc_all, gc_exon_all, at_exon_all],
                         ["otherGC-exons", "GC-exons", "AT-exons"], output)


def control_query_execution(cnx, exon_type):
    """
    query launcher

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return: (list of tuple)
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT t1.gene_id, t1.exon_pos,
                              t1.iupac_exon, t1.upstream_intron_size, 
                              t1.downstream_intron_size
                       FROM sed t1
                       WHERE t1.exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT t1.gene_id, t1.exon_pos,
                              t1.iupac_exon, t1.upstream_intron_size, 
                              t1.downstream_intron_size
                       FROM sed t1
                    """
    cursor.execute(query)
    return cursor.fetchall()


def get_control_exon_information(cnx, exon_type, exon2remove):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` \
    exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of list of 2int) list of exon to remove from\
     the control list of exons
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    result = control_query_execution(cnx, exon_type)
    # turn tuple into list
    print("number of control exon before removing bad ones : %s" % len(result))
    nresult = [list(exon) for exon in result if
               [exon[0], exon[1]] not in exon2remove]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    min_flanking_intron_size = np.median([np.nanmin(
                                    np.array([exon[3], exon[4]], dtype=float))
                                          for exon in nresult])
    gc_content = np.median([float(exon[2].split(";")[4]) for exon in
                           nresult])
    print("median GC content : %s" % gc_content)
    print("median min flaking intron size : %s" % min_flanking_intron_size)
    return min_flanking_intron_size, gc_content


def get_exons_of_interest(cnx, exon_type, exon2remove, min_intron_size,
                          gc_content, mtype):
    """

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon2remove: (list of list of 2int) list of exon to remove from\
     the control list of exons
    :param min_intron_size: (int) the median size of the smallest flaking \
    intron
    :param gc_content: (int) the median value of gc content
    :param mtype: gc or at
    :return: (list of 2 int) list of exons
    """
    result = control_query_execution(cnx, exon_type)
    print("number of control exon before removing bad ones : %s" % len(result))
    if mtype == "GC":
        nresult = [[exon[0], exon[1]] for exon in result if
                   [exon[0], exon[1]] not in exon2remove and
                   np.nanmin(np.array([exon[3], exon[4]], dtype=float)) <
                   min_intron_size and
                   float(exon[2].split(";")[4]) > gc_content]
    else:
        nresult = [[exon[0], exon[1]] for exon in result if
                   [exon[0], exon[1]] not in exon2remove and
                   np.nanmin(np.array([exon[3], exon[4]], dtype=float)) >
                   min_intron_size and
                   float(exon[2].split(";")[4]) < gc_content]
    print("number of control exon after removing bad ones : %s" % len(nresult))
    return nresult


def find_first_exons(cnx, exon_type):
    """
    Find every first exons with the type ``exon_type``.

    :param cnx: (sqlite3 connection object) connection to sed database
    :param exon_type: (srt) the control exon type chosen (CCE, ACE)..
    :return: (list of 2 int) list of first exons
    """
    cursor = cnx.cursor()
    query = """SELECT t1.gene_id, t1.exon_pos
               FROM sed t1
               WHERE t1.exon_type LIKE '%{}%'
               AND t1.exon_pos = 1""".format(exon_type)
    cursor.execute(query)
    res = cursor.fetchall()
    res = [list(exon) for exon in res]
    print('First exons :', res[0:5])
    return


def find_gene_id(cnx, exon_type):
    """
    Find every gene id.

    :param cnx: (sqlite3 connection object) connection to sed database
    :param exon_type: (srt) the control exon type chosen (CCE, ACE)..
    :return: (list of int) list genes
    """
    cursor = cnx.cursor()
    query = """SELECT t1.gene_id
               FROM sed t1
               WHERE t1.exon_type LIKE '%{}%'""".format(exon_type)
    cursor.execute(query)
    res = cursor.fetchall()
    return [exon[0] for exon in res]


def find_last_exons(cnx, exon_type):
    """
    Find every last intron.

    :param cnx: (sqlite3 connection object) connection to sed database
    :param exon_type: (srt) the control exon type chosen (CCE, ACE)..
    :return:
    """
    gene_list = find_gene_id(cnx, exon_type)
    result = []
    cursor = cnx.cursor()
    for mgene in gene_list:
        query = """SELECT t1.gene_id, t1.exon_pos
               FROM sed t1
               WHERE t1.gene_id = {}
               ORDER BY t1.exon_pos DESC
               LIMIT 1""".format(mgene)
        cursor.execute(query)
        res = cursor.fetchone()
        exon = "_".join(list(map(str, res)))
        result.append(exon)
    result = list(np.unique(result))
    print("Last exons :", result[0:5])
    return result


def create_unregulated_exon_list(cnx, output, exon_type):
    """
    Create the list of GC/At rich unregulated exons

    :param cnx: (sqlite3 connect objecy) connection to sed database
    :param output: (str) path were the exon files will be created
    :param exon_type: (str) the type of control exons
    """
    gc_exon_all = get_exons_list(cnx, group_factor.gc_rich_down, "down")
    at_exon_all = get_exons_list(cnx, group_factor.at_rich_down, "down")
    gc_exon = [exon for exon in gc_exon_all if exon not in at_exon_all]
    at_exon = [exon for exon in at_exon_all if exon not in gc_exon_all]
    first_exons = find_first_exons(cnx, exon_type)
    last_exons = find_last_exons(cnx, exon_type)
    exon2remove = gc_exon + at_exon + last_exons
    print("exon to remove : %s" % len(exon2remove))
    min_intron_size, gc_content = \
        get_control_exon_information(cnx, exon_type, exon2remove)
    gc_exon = get_exons_of_interest(cnx, exon_type, exon2remove,
                                    min_intron_size, gc_content, "GC")
    at_exon = get_exons_of_interest(cnx, exon_type, exon2remove,
                                    min_intron_size, gc_content, "AT")
    gc_exon = [exon for exon in gc_exon if exon not in at_exon]
    at_exon = [exon for exon in at_exon if exon not in gc_exon]
    print("gc_exon : %s exons" % len(gc_exon))
    print("at_exon : %s exons" % len(at_exon))
    file_writer(gc_exon, "GC_unregulated", output)
    file_writer(at_exon, "AT_unregulated", output)
    cnx.close()


def main():
    exon_type = "CCE"
    seddb = os.path.realpath(os.path.dirname(os.path.dirname(__file__)
                                             ).replace("src","data/sed.db"))
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(
        os.path.dirname(os.path.dirname(__file__)
                        )).replace("src", "result/figure_2_3/exon_list/")

    if not os.path.isdir(output):
        os.mkdir(output)
    # print("%sCreate GA and CT exons list%s" % ("-" * 20, "-" * 20))
    # create_ct_ga_rich_exon_list(cnx, output)
    # print("%sCreate other GC exons list%s" % ("-" * 20, "-" * 20))
    # create_othergc_exon_file(cnx, output)
    print("%sCreate unregulated AT/GC exons list%s" % ("-" * 20, "-" * 20))
    create_unregulated_exon_list(cnx, output, exon_type)
    cnx.close()
    print("Finished !")


if __name__ == "__main__":
    main()
