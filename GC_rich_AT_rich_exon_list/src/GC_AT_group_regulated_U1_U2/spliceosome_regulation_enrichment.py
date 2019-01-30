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
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__).replace("/GC_AT_group_regulated_U1_U2", "")))
import group_factor
from figure_creator import get_exons_list


def exon_intersection(exon_list1, exon_list2):
    """
    Return the intersection of 2 exon lists

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the intersection between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon in exon_list2]
    return exon_list


def exon_difference(exon_list1, exon_list2):
    """
    Return the difference of 2 exon lists

    :param exon_list1: (list of 2 int) gene_id and exon_position
    :param exon_list2: (list of 2 int) gene_id and exon position
    :return: (list of 2 int) gene_id and exon position, the difference between ``exon_list1`` and ``exon_list2``
    """
    # intersection
    exon_list = [exon for exon in exon_list1 if exon not in exon_list2]
    return exon_list


def frequency_test(obs1, tot1, obs2, tot2):
    """
    Chiq test

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
    Get the wanted control exons
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
    Make a subsampling of the data in ``my_list``

    :param my_list: (list of list of 2 int) list of exons
    :param nb: (int) the number of value we want to subsample in ``my_list``
    :return: (list of list of 2 int) ``my_list`` subsampled
    """
    list_values = np.random.choice(range(len(my_list)), nb, replace=False)
    subsampled_list = []
    for val in list_values:
        subsampled_list.append(my_list[val])
    return subsampled_list


def analysis_maker(gc_exon_list, at_exon_list, gc_pure_exon_list, at_pure_exon_list, splicesome_dic,
                   u1_u2_intersection, nb_iteration):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list, one containing exons regulated by splicing factors regulating AT rich exons and the other \
    containing exons regulated by splicing factors regulating GC rich exons
    :param gc_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating GC rich exons.
    :param at_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating AT rich exons.
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
    super_dict = {}
    at_gc = ["pure"]
    ctrl = ["all"]
    tot = "NA"
    for gc_at_set in at_gc:
        if gc_at_set == "all":
            cur_gc_list = gc_exon_list
            cur_at_list = at_exon_list
        else:
            cur_gc_list = gc_pure_exon_list
            cur_at_list = at_pure_exon_list
        for ctrl_set in ctrl:
            for key in splicesome_dic.keys():
                if ctrl_set == "all":
                    cur_ctrl = splicesome_dic[key]
                else:
                    cur_ctrl = exon_difference(splicesome_dic[key], u1_u2_intersection)
                list_pval = []
                list_gc_intersect = []
                list_at_intersect = []
                for i in range(nb_iteration):
                    sys.stdout.write("%s/%s     \r" % (i, nb_iteration))
                    # print(i, len(cur_gc_list), len(cur_at_list))
                    if len(cur_gc_list) > len(cur_at_list):
                        new_gc_list = subsample(cur_gc_list, len(cur_at_list))
                        new_at_list = cur_at_list
                    else:
                        new_at_list = subsample(cur_at_list, len(cur_gc_list))
                        new_gc_list = cur_gc_list
                    gc_intersect = len(exon_intersection(new_gc_list, cur_ctrl))
                    list_gc_intersect.append(gc_intersect)
                    at_intersect = len(exon_intersection(new_at_list, cur_ctrl))
                    list_at_intersect.append(at_intersect)
                    tot = len(cur_ctrl)
                    pval = frequency_test(gc_intersect, tot, at_intersect, tot)
                    list_pval.append(pval)

                analysis_dic["GC-AT-%s-vs-%s-%s" % (gc_at_set, key, ctrl_set)] = \
                    [np.mean(list_pval), np.mean(list_gc_intersect), tot, np.mean(list_at_intersect),
                     tot, float(sum(np.array(list_pval) > 0.05) / len(list_pval))]
                super_dict["GC-AT-%s-vs-%s-%s" % (gc_at_set, key, ctrl_set)] = \
                    [list_pval, list_gc_intersect, tot, list_at_intersect, tot]
    return analysis_dic, super_dict


def analysis_maker_ctrl(gc_exon_list, gc_pure_exon_list, ctrl_exons, splicesome_dic,
                        u1_u2_intersection, nb_iteration, exon_type, test_set):
    """
    Make the comparison analysis of the frequencies of exons regulated by every spliceosome(U1) factors from \
    two exons list, one containing exons regulated by splicing factors regulating AT rich exons and the other \
    containing exons regulated by splicing factors regulating GC rich exons
    :param gc_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating GC rich exons.
    :param gc_pure_exon_list: (list of list of 2 int) list of exons regulated by splicing factors regulating \
    GC rich exons and not regulated by splicing factors regulating AT rich exons.
    :param ctrl_exons: (list of list of 2int) list of control exons
    :param splicesome_dic: (dict of list of 2 int) dictionary that contains every list of exons regulated by a \
    particular U1 spliceosome factors.
    :param u1_u2_intersection: (list of list of 2 int) list of exons regulated by both u1 and u2 factors.
    :param nb_iteration: (int) the number of iteration we gonna make
    :param exon_type: (string) the type of exon studied
    :param test_set: (string) test set
    :return:
        - analysis_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done to \
    compare two different list of exons.
        - super_dic: (dictionary of list of int) the dictionary that contains the enrichment analysis done for every \
        subsampling realized.
    """
    analysis_dic = {}
    super_dict = {}
    at_gc = ["pure"]
    ctrl = ["all"]
    tot = "NA"
    for gc_at_set in at_gc:
        if gc_at_set == "all":
            cur_gc_list = gc_exon_list
        else:
            cur_gc_list = gc_pure_exon_list
        for ctrl_set in ctrl:
            for key in splicesome_dic.keys():
                if ctrl_set == "all":
                    cur_ctrl = splicesome_dic[key]
                else:
                    cur_ctrl = exon_difference(splicesome_dic[key], u1_u2_intersection)
                list_pval = []
                list_gc_intersect = []
                list_ctrl_intersect = []
                for i in range(nb_iteration):
                    sys.stdout.write("%s/%s     \r" % (i, nb_iteration))
                    # print(i, len(cur_gc_list), len(cur_at_list))
                    new_ctrl = subsample(ctrl_exons, len(cur_gc_list))
                    new_gc_list = cur_gc_list
                    gc_intersect = len(exon_intersection(new_gc_list, cur_ctrl))
                    list_gc_intersect.append(gc_intersect)
                    ctrl_intersect = len(exon_intersection(new_ctrl, cur_ctrl))
                    list_ctrl_intersect.append(ctrl_intersect)
                    tot = len(cur_ctrl)
                    pval = frequency_test(gc_intersect, tot, ctrl_intersect, tot)
                    list_pval.append(pval)

                analysis_dic["%s-%s-%s-vs-%s-%s" % (test_set, exon_type, gc_at_set, key, ctrl_set)] = \
                    [np.mean(list_pval), np.mean(list_gc_intersect), tot, np.mean(list_ctrl_intersect),
                     tot, float(sum(np.array(list_pval) > 0.05) / len(list_pval))]
                super_dict["%s-%s-%s-vs-%s-%s" % (test_set, exon_type, gc_at_set, key, ctrl_set)] = \
                    [list_pval, list_gc_intersect, tot, list_ctrl_intersect, tot]
    return analysis_dic, super_dict


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
    with open("%sresult_enrichment_analysis_%s_iteration_%s_vs_%s.txt" %
              (output, nb_iteration, test_set, comparison), "w") as outfile:
        outfile.write("# nb iteration : %s" % nb_iteration)
        outfile.write("Name_analysis\tratio1\tcalcul_ratio1\tratio2\tcalcul_ratio2"
                      "\tpval_on_mean_count\tmean_pval\tfreq_pval>0.05\n")
        for key in sorted(analysis_dic.keys()):
            name_gc = "%s(%s)" % (key.split("-")[0], key.split("-")[2])
            name_at = "%s(%s)" % (key.split("-")[1], key.split("-")[2])
            name_ctrl = "(".join(key.split("-")[4:]) + ")"
            pval = frequency_test(round(analysis_dic[key][1]), int(analysis_dic[key][2]),
                                  round(analysis_dic[key][3]), int(analysis_dic[key][4]))
            outfile.write("%s\t%s/%s\t(%s-%s/%s)\t%s/%s\t(%s-%s/%s)\t%s\t%s\t%s\n" %
                          (key, analysis_dic[key][1], analysis_dic[key][2], name_gc,
                           name_ctrl, name_ctrl, analysis_dic[key][3], analysis_dic[key][4], name_at,
                           name_ctrl, name_ctrl, pval, analysis_dic[key][0], analysis_dic[key][5]))


def main(test_set, comparison):
    """
    Make the enrichment analysis comparing the frequencies of exon regulated by splicesome factors \
     for an AT and GC exons list.

    :param comparison: (string) the comparison set of exon used
    :param test_set: (string) the test set used
    """
    nb_iteration = 10
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/GC_AT_group_regulated_U1_U2",
                                                                 "result/GC_AT_group_regulated_U1_U2/")
    if not os.path.isdir(output):
        os.mkdir(output)
    # if test_set == "GC":
    #     div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
    #                  "U1" : group_factor.u1_factors, "U1_DDX5": list(group_factor.u1_factors) + ["DDX5_DDX17"],
    #                  "SNRPC": ["SNRPC"], "SNRNP70": ["SNRNP70"], "DDX5_17": ["DDX5_DDX17"]}
    # else:
    #     div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
    #                  "U2" : group_factor.u2_factors, "SF1": ["SF1"], "U2AF1": ["U2AF1"], "U2AF2": ['U2AF2'],
    #                  "SF3B1": ["SF3B1"], "SF3B4": ["SF3B4"]}
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
    if comparison == "AT" or comparison == "GC":
        analysis_dic, super_dict = analysis_maker(dic_exon["GC_rich"], dic_exon["AT_rich"],
                                                  dic_exon["GC_pure"], dic_exon["AT_pure"],
                                                  dic_spliceosome, u1_u2_intersection, nb_iteration)
    else:
        ctrl_exons = get_control_exon(cnx, comparison)
        if test_set == "GC":
            analysis_dic, super_dict = analysis_maker_ctrl(dic_exon["GC_rich"], dic_exon["GC_pure"],
                                                           ctrl_exons, dic_spliceosome, u1_u2_intersection,
                                                           nb_iteration, comparison, test_set)
        else:
            analysis_dic, super_dict = analysis_maker_ctrl(dic_exon["AT_rich"], dic_exon["AT_pure"], ctrl_exons,
                                                           dic_spliceosome, u1_u2_intersection, nb_iteration,
                                                           comparison, test_set)
    file_writer(output, analysis_dic, nb_iteration, test_set, comparison)
    with open("%ssuper_dict_%s_%s.txt" % (output, nb_iteration, comparison), "w") as outfile:
        outfile.write("super_dict=%s\n" % str(super_dict))
    cnx.close()


if __name__ == "__main__":
    if sys.argv[1] in ["GC", "AT"] and sys.argv[2] in ["GC", "AT", "CCE"]:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Wrong parameter values")
