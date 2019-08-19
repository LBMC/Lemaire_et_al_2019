#!/usr/bin/env python3.5

"""
Description:
    The goal of this script is to produce a figure that show the length of \
    control fasterdb exons and compare the variance in GC content of short \
    and long exons.
"""

import sqlite3
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import union_dataset_function as udf
from scipy.stats import levene


nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7, "K": 8, "M": 9}
dnt_dic = {"AA": 0, "AC": 1, "AG": 2, "AT": 3, "CA": 4, "CC": 5,
           "CG": 6, "CT": 7, "GA": 8, "GC": 9, "GG": 10, "GT": 11,
           "TA": 12, "TC": 13, "TG": 14, "TT": 15}


def get_list_of_value(cnx, exon_list, target_column):
    """
    Get the individual values for ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds \
    to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get \
    information on exons.
    :return: (dict of float) values of ``target_column`` for the exons in \
     ``exon_list``.
    """
    gene_fields = ["gene_size", "nb_intron_gene", "median_intron_size",
                   "iupac_gene", "dnt_gene"]
    cursor = cnx.cursor()
    res = {"exon": [], target_column: []}
    if target_column not in gene_fields:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % \
                    (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                res["exon"].append(exon)
                res[target_column].append(r)
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                           FROM sed
                           where gene_id = %s
                           AND exon_pos = %s """ % \
                        (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                res["exon"].append(exon)
                res[target_column].append(r)
                redundancy_gene_dic[exon[0]] = 1
    return res


def get_control_exon(cnx, exon_type, exon_2_remove, regulation):
    """
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :param exon_2_remove: (string) exons regulated by a splicing factor
    :param regulation: (string)  up or down
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos
                   FROM sed
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    print("%s exons : %s" % (exon_type, len(result)))
    nresult = [list(exon) for exon in result
               if [exon[0], exon[1]] not in exon_2_remove]
    print("%s exons not %s-regulated by splicing factors : %s" %
          (exon_type, regulation, len(nresult)))
    return nresult


def make_histogram(list_value, output, figname, scale, log=False):
    """
    From a list of values, create an histogram.

    :param list_value: (list of float/int) list of values
    :param outout: (str) location were the figure will be created
    :param figname: (str) the name of the figure
    :param scale: (str) the scale chosen
    :param log: (boolean) True to have an x-axis in log scale, false else
    """
    sns.set()
    sns.set_context("poster")
    fig = plt.figure(figsize=(10, 10))
    if log:
        bins = np.logspace(1, 3, 20)
    else:
        bins = 20
    ax = sns.distplot(list_value, axlabel=scale.replace("_", " "), hist_kws=dict(alpha=1),
                      rug=False, kde=False, bins=bins)
    if log:
        fig.get_axes()[0].set_xscale("log")
    ax.set_title("%s distribution" % scale.replace("_", " "))
    plt.savefig("%s/%s.pdf" % (output, figname))
    plt.clf()
    plt.cla()
    plt.close()


def get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt_dnt: (string) a nucleotide or di_nucleotide
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                if len(nt_dnt) == 1:
                    res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                else:
                    res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                if r is not None:
                    if len(nt_dnt) == 1:
                        res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                    else:
                        res.append(float(r.split(";")[dnt_dic[nt_dnt]]))
                redundancy_gene_dic[exon[0]] = 1
    return res


def get_two_groups_of_exon(dic_size, threshold, target_columns):
    """
    Return two list of exons, one having size below or equal to ``threshold`` \
    and the other having a size greater than ``threshold``.

    :param dic_size: (dictionary of list of float) dictionary of exons size
    :param threshold: (int) the size chosen to create the list of exons
    :param target_columns: (str) exon_size
    :return: (2 lists of of list of 2 int) list of exons size
    """
    small_exons = []
    big_exons = []
    for i in range(len(dic_size["exon"])):
        if dic_size[target_columns][i] <= threshold:
            small_exons.append(dic_size["exon"][i])
        else:
            big_exons.append(dic_size["exon"][i])
    return small_exons, big_exons


def write_test_result(big_gc, small_gc, threshold, output, name_exons):
    """
    Write a file containing the result of a levene test between the GC-content \
     of a group of big exons and a group of small exons.

    :param big_gc: (list of float) list of gc content of big exons : exons \
    with a size greater than ``threshold``
    :param small_gc: (list of float) list of gc content of big exons : exons \
    with a size lower or equal than ``threshold``
    :param threshold: (float) the size chosen to separate big and small exons
    :param output: (str) folder where the result will be created
    :param name_exons: (str ) the name of the exons studied
    """
    w, pvalue = levene(small_gc, big_gc)
    print("levene (homoscedaticity test) pvalue = %s" % pvalue)
    with open("%s/pvalue_levene_test_%s.txt" % (output, name_exons), "w") \
            as ouf:
        ouf.write("Levene test comparing exons:\n")
        ouf.write("\t-with a size below/equal to %s nt (%s exons)\n" %
                  (threshold, len(small_gc)))
        ouf.write("\t-with a size greater to %s nt (%s exons)\n" %
                  (threshold, len(big_gc)))
        ouf.write("====  PVALUE : %s =====\n" % pvalue)
        ouf.write("std gc content small exons : %s\n" % (np.std(small_gc)))
        ouf.write("std gc content big exons : %s\n" % (np.std(big_gc)))


def my_level_analysis(cnx, exon_type, output, regulation, size_threshold,
                      target_column, level="control"):
    """
    Create the histogram of the size of exons and make a levene test \
    to test if the variance of GC content of a group of big exons \
    and a group of small exons is different.

    :param cnx: (pymysql connection object) connection to Sed database.
    :param exon_type: (str) the type of control exons to analyse
    :param output: (str) folder where the results will be created
    :param regulation: (str) the regulation
    :param size_threshold: (int) the threshold
    :param target_column: (str) the feature of interest
    :param level: (str) the level
    """
    sizefig = "hist_of_%s_exon_size" % exon_type
    if level == "control":
        exon_2_remove = udf.get_exon_regulated_by_sf(cnx, regulation)
        exon_list = get_control_exon(cnx, exon_type, exon_2_remove, regulation)
    else:
        exon_list = udf.get_exon_regulated_by_sf(cnx, regulation)
        exon_type = "SF-down"
    dic_size = get_list_of_value(cnx, exon_list, target_column)
    # dic_size = {"exon": [], target_column: []}
    # for i in range(len(tmp["exon"])):
    #     if tmp[target_column][i] > 10:
    #         dic_size["exon"].append(tmp["exon"][i])
    #         dic_size[target_column].append(tmp[target_column][i])
    list_size = np.array(dic_size["exon_size"])
    print(" min : %s" % min(list_size))
    print(" max : %s" % max(list_size))
    print(" nb exons having a size below/equal %s : %s" %
          (size_threshold, len(list_size[list_size <= size_threshold])))
    make_histogram(list_size, output, sizefig, target_column, log=True)
    small_exons, big_exons = get_two_groups_of_exon(dic_size, size_threshold,
                                                    target_column)
    print(" nb exons having a size below/equal to %s nt : %s" %
          (size_threshold, len(small_exons)))
    print(" nb exons having a size greater to %s nt : %s" %
          (size_threshold, len(big_exons)))
    small_gc = get_list_of_value_iupac_dnt(cnx, small_exons,
                                           "iupac_exon", "S")
    big_gc = get_list_of_value_iupac_dnt(cnx, big_exons,
                                           "iupac_exon", "S")
    make_histogram(small_gc, output, "gc_content_small_%s_exon" % exon_type,
                   "GC content exons <= %s nt" % size_threshold)
    make_histogram(big_gc, output, "gc_content_big_%s_exon" % exon_type,
                   "GC content exons > %s nt" % size_threshold)
    write_test_result(big_gc, small_gc, size_threshold, output, exon_type)


def main():
    """
    Create the graphics wanted
    """
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    seddb = base + "/data/sed.db"
    output = base + "/result/variance_analysis"
    if not os.path.isdir(output):
        os.mkdir(output)
    exon_type = "CCE"
    regulation = "down"
    target_column = "exon_size"
    size_threshold = 27
    cnx = sqlite3.connect(seddb)
    my_level_analysis(cnx, exon_type, output, regulation, size_threshold,
                           target_column)
    my_level_analysis(cnx, exon_type, output, regulation, size_threshold,
                           target_column, level="SF-down")
    cnx.close()


if __name__ == "__main__":
    main()