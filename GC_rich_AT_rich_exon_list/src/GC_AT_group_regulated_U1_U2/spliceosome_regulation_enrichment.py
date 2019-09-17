#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:

    This program will see if a list of exons GC rich is more often regulated by u1 than AT-rich or CCE:ACE... list\
    of exons.
"""

import sys
import numpy as np
import os
import sqlite3
import rpy2.robjects as robj
import rpy2.robjects.vectors as v
from barplot_pvalue_maker import fig_3g
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__).replace("/GC_AT_group_regulated_U1_U2", "")))
import group_factor
from figure_creator import get_exons_list


def exon_intersection(exon_list1, exon_list2):
    """
    Return the intersection of 2 exon lists.

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the intersection between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon in exon_list2]
    return exon_list


def exon_difference(exon_list1, exon_list2):
    """
    Return the difference of 2 exon lists.

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the difference between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon not in exon_list2]
    return exon_list


def frequency_test(obs1, tot1, obs2, tot2):
    """
    Chiq test.

    :param obs1: (int) the count number of an amino acid X in the set of protein 1.
    :param tot1: (int) the total number of amino acids in the set of protein 1.
    :param obs2: (int) the count number of an amino acid X in the set of protein 2.
    :param tot2: (int) the total number of amino acids in the set of protein 2.
    :return: proportion test p-value
    """

    chisq = robj.r("""

        function(vect){
            m<-matrix(vect, byrow=T, nrow=2)
            return(chisq.test(m)$p.value)
        }

                       """)
    rm1 = tot1 - obs1
    rm2 = tot2 - obs2
    vect = v.FloatVector([obs1, rm1, obs2, rm2])
    pval = float(chisq(vect)[0])
    if np.isnan(pval):
        return "NA"
    else:
        return pval


def get_control_exon(cnx, exon_type):
    """
    Get the wanted control exons.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return:
        * result: (list of 2 int) list of control exons
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos
                   FROM sed"""
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    values = []
    for val in result:
        values.append(list(map(str, list(val))))
    return values


def subsample(my_list, nb):
    """
    Make a subsampling of the data in ``my_list``.

    :param my_list: (list of list of 2 int) list of exons
    :param nb: (int) the number of value we want to subsample in ``my_list``
    :return: (list of list of 2 int) ``my_list`` subsampled
    """
    list_values = np.random.choice(range(len(my_list)), nb, replace=False)
    subsampled_list = []
    for val in list_values:
        subsampled_list.append(my_list[val])
    return subsampled_list


def calculate_index(value, list_of_values):
    """
    :param value: a value corresponding to a codon frequency
    :param list_of_values: a list of codon frequencies
    :return: an empirical "p-value" calculated by the position of the value "value" in the list_of_value \
    sorted in the ascending order
    """
    list_of_values.sort()
    val1 = 0
    val2 = 1
    for i in range(len(list_of_values)):  # proportion of values in the list_of_values bigger than the value "value"
        if value <= list_of_values[i]:
            val1 = float(len(list_of_values) - i) / len(list_of_values)
            break  # as soon as we find the index, we quit the loop
    for i in range(len(list_of_values)):  # proportion of values in the list_of_values smaller than the value "value"
        if value < list_of_values[i]:
            val2 = float(i) / len(list_of_values)
            break  # as soon as we find the index, we quit the loop*
    if val2 == -1 and list_of_values[-1] == value:
        return 1, "="
    if val1 != -1 and val2 != -1:
        if val1 < val2:
            return min(val1, 1), "+"
        if val1 == val2:
            return min(val1, 1), "="
        else:
            return min(val2, 1), "-"
    elif val1 == -1 and val2 == -1:
        return "NA", "NA"
    elif val1 == -1 and val2 != -1:
        return 0, "+"
    else:
        return 0, "-"


def adapt_reg(pvalue, regulation, alpha=0.05):
    """
    Return the adapted regulation.

    :param pvalue: (float) the p-value
    :param regulation: (string) + for enrichment - for impoverishement = for nothing
    :param alpha: (string)
    :return: (string) return + or - if the pvalue is significant
    """
    if pvalue <= alpha:
        return regulation
    return "="


def analysis_maker(gc_pure_exon_list, at_pure_exon_list, splicesome_dic,
                   u1_u2_intersection, nb_iteration):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list, one containing exons regulated by splicing factors regulating AT rich exons and the other \
    containing exons regulated by splicing factors regulating GC rich exons.

    :param gc_pure_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating \
    GC rich exons and not regulated by splicing factors regulating AT rich exons.
    :param at_pure_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating \
    AT rich exons and not regulated by splicing factors regulating GC rich exons.
    :param splicesome_dic: (dict of list of 2 int) dictionary that contains every list of exons regulated by a \
    particular U1 spliceosome factors.
    :param u1_u2_intersection: (list of list of 2 int) list of exons regulated by both u1 and u2 factors.
    :param nb_iteration: (int) the number of iteration we gonna make
    :return:
        - analysis_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done to \
    compare two different list of exons.
        - super_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done for every \
        subsampling realized.
    """
    analysis_dic = {}
    super_dic = {}
    ctrl_set = "all"

    cur_gc_list = gc_pure_exon_list
    cur_at_list = at_pure_exon_list
    for key in splicesome_dic.keys():
        if ctrl_set == "all":
            cur_ctrl = splicesome_dic[key]
        else:
            cur_ctrl = exon_difference(splicesome_dic[key], u1_u2_intersection)
        list_at_intersect = []
        gc_intersect = len(exon_intersection(cur_gc_list, cur_ctrl))
        for i in range(nb_iteration):
            sys.stdout.write("%s/%s     \r" % (i, nb_iteration))
            new_at_list = None
            if len(cur_gc_list) > len(cur_at_list):
                print("Error : the GC-exons should not be greater than AT-exons list !")
                exit(1)
            else:
                new_at_list = subsample(cur_at_list, len(cur_gc_list))
            at_intersect = len(exon_intersection(new_at_list, cur_ctrl))
            list_at_intersect.append(at_intersect)
        pval, reg = calculate_index(gc_intersect, list_at_intersect)
        pval = max(pval, float(1 / nb_iteration))
        analysis_dic["GC-AT-vs-%s-%s" % (key, ctrl_set)] = \
            [pval, adapt_reg(pval, reg)]
        super_dic["GC-AT-vs-%s-%s" % (key, ctrl_set)] = \
            [pval, gc_intersect, list_at_intersect]
    return analysis_dic, super_dic


def get_exon_list_file(list_file):
    """

    :param list_file: (str) a file containing exons
    :return: (list of 2 int) a list of exons
    """
    exon_list = []
    with open(list_file, "r") as infile:
        for line in infile:
            line = line.replace("\n", "")
            line = line.split("\t")
            exon_list.append(line)
    return exon_list


def analysis_maker_bis(list_file, name_file, splicesome_dic, nb_iteration):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list.

    :param list_file: (list of str) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (list of str) the name of each files of exons \
    given in ``list_file``
    :param splicesome_dic: (dict of list of 2 int) dictionary that contains every list of exons regulated by a \
    particular U1 spliceosome factors.
    :param u1_u2_intersection: (list of list of 2 int) list of exons regulated by both u1 and u2 factors.
    :param nb_iteration: (int) the number of iteration we gonna make
    :return:
        - analysis_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done to \
    compare two different list of exons.
        - super_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done for every \
        subsampling realized.
    """
    analysis_dic = {}
    super_dic = {}
    if len(list_file) != 2:
        raise IndexError("list_file should have a length of 2 elements")
    cur_list1 = get_exon_list_file(list_file[0])
    name_list1 = name_file[0].replace("-", "_")
    print("%s : %s" % (name_list1, len(cur_list1)))
    cur_list2 = get_exon_list_file(list_file[1])
    name_list2 = name_file[1].replace("-", "_")
    print("%s : %s" % (name_list2, len(cur_list2)))
    for key in splicesome_dic.keys():
        cur_ctrl = splicesome_dic[key]
        list_intersect = []
        if len(cur_list1) < len(cur_list2):
            intersect1 = len(exon_intersection(cur_list1, cur_ctrl))
            interb = len(exon_intersection(cur_list2, cur_ctrl))
        else:
            intersect1 = len(exon_intersection(cur_list2, cur_ctrl))
            interb = len(exon_intersection(cur_list1, cur_ctrl))
        for i in range(nb_iteration):
            sys.stdout.write("%s/%s     \r" % (i, nb_iteration))
            if len(cur_list1) > len(cur_list2):
                intersect2 = subsample(cur_list1, len(cur_list2))
            else:
                intersect2 = subsample(cur_list2, len(cur_list1))
            at_intersect = len(exon_intersection(intersect2, cur_ctrl))
            list_intersect.append(at_intersect)
        pval, reg = calculate_index(intersect1, list_intersect)
        pval = max(pval, float(1 / nb_iteration))
        name_file = [name_list1, name_list2]
        if len(cur_list2) < len(cur_list1):
            name_file = name_file[::-1]

        analysis_dic["%s-%s-vs-%s" % (name_file[0], name_file[1], key)] = \
            [pval, reg]
        super_dic["%s-%s-vs-%s" % (name_file[0], name_file[1], key)] = \
            [pval, intersect1, list_intersect, interb]
    return analysis_dic, super_dic


def file_writer(output, analysis_dic, nb_iteration, test_set, comparison):
    """
    Write the enrichment analysis file.

    :param output: (string) path were the result will be  created
    :param analysis_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done to \
    compare two different list of exons.
    :param nb_iteration: (int) the number of iteration used to get the ``analysis_dic``
    :param comparison: (string) list of exons compared to the gc list of exons
    :param test_set: (string) the test set used
    """
    filename = "%sresult_enrichment_analysis_%s_iteration_%s_vs_%s.txt" % \
              (output, nb_iteration, test_set, comparison)
    with open(filename, "w") as outfile:
        outfile.write("# nb iteration : %s" % nb_iteration)
        outfile.write("Name_analysis\tenpirical-pvalue\tregulation\n")
        for key in sorted(analysis_dic.keys()):
            outfile.write("%s\t%s\t%s\n" % (key, analysis_dic[key][0], analysis_dic[key][1]))
    return filename


def main():
    """
    Make the enrichment analysis comparing the frequencies of exon regulated by splicesome factors \
     for an AT and GC exons list.
    """
    nb_iteration = 10
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2",
                                                                 "result/GC_AT_group_regulated_U1_U2/")
    if not os.path.isdir(output):
        os.mkdir(output)
    div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
                 "SNRPC": ["SNRPC"], "SNRNP70": ["SNRNP70"], "DDX5_17": ["DDX5_DDX17"],
                 "SF1": ["SF1"], "U2AF1": ["U2AF1"], "U2AF2": ['U2AF2'], "SF3A3": ["SF3A3"], "SF3B4": ["SF3B4"]}
    dic_exon = {}
    for name_group in div_group.keys():
        print("Getting all exon regulated by %s factor" % name_group)
        dic_exon[name_group] = get_exons_list(cnx, div_group[name_group], "down")
    at_gc_intersection = exon_intersection(dic_exon["AT_rich"], dic_exon["GC_rich"])
    u1_u2_intersection = exon_intersection(get_exons_list(cnx, group_factor.u1_factors, "down"),
                                           get_exons_list(cnx, group_factor.u2_factors, "down"))
    print("GC-AT group intersection : %s exons" % len(at_gc_intersection))
    print("U1-U2 interesection : %s exons" % len(u1_u2_intersection))
    dic_exon["GC_pure"] = exon_difference(dic_exon["GC_rich"], at_gc_intersection)
    dic_exon["AT_pure"] = exon_difference(dic_exon["AT_rich"], at_gc_intersection)
    for key in dic_exon:
        print("%s : %s" % (key, len(dic_exon[key])))
    dic_spliceosome = {}
    for key in dic_exon.keys():
        if "AT" not in key and "GC" not in key:
            dic_spliceosome[key] = dic_exon[key]

    analysis_dic, super_dict = analysis_maker(dic_exon["GC_pure"], dic_exon["AT_pure"],
                                              dic_spliceosome, u1_u2_intersection, nb_iteration)
    file_writer(output, analysis_dic, nb_iteration, "GC", "AT")
    with open("%ssuper_dict_%s.txt" % (output, nb_iteration), "w") as outfile:
        outfile.write("super_dict=%s\n" % str(super_dict))
    cnx.close()


def main_3g(list_file, name_file, seddb, output, reverse):
    """
    Make the enrichment analysis comparing the frequencies of exon regulated by splicesome factors \
     for an AT and GC exons list.
    """
    nb_iteration = 10000
    cnx = sqlite3.connect(seddb)
    div_group = {"SNRPC": ["SNRPC"], "SNRNP70": ["SNRNP70"],
                 # "DDX5_17": ["DDX5_DDX17"],
                 "SF1": ["SF1"], "U2AF1": ["U2AF1"], "U2AF2": ['U2AF2'],
                 # "SF3A3": ["SF3A3"], "SF3B4": ["SF3B4"]
                 }
    dic_exon = {}
    for name_group in div_group.keys():
        print("Getting all exon regulated by %s factor" % name_group)
        dic_exon[name_group] = get_exons_list(cnx, div_group[name_group], "down")
    for key in dic_exon:
        print("%s : %s" % (key, len(dic_exon[key])))
    dic_spliceosome = {}
    for key in dic_exon.keys():
        if "AT" not in key and "GC" not in key:
            dic_spliceosome[key] = dic_exon[key]

    analysis_dic, super_dict = analysis_maker_bis(list_file, name_file,
                                                  dic_spliceosome, nb_iteration)
    filename = file_writer(output, analysis_dic, nb_iteration, name_file[0],
                name_file[1])
    fig_3g(filename, output, reverse)
    with open("%ssuper_dict_%s.txt" % (output, nb_iteration), "w") as outfile:
        outfile.write("super_dict=%s\n" % str(super_dict))
    cnx.close()

if __name__ == "__main__":
    main()
