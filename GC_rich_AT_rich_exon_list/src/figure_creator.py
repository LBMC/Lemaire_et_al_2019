#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:

    The goal of this script is to create files containing the list of exon that we will studied and creates venn \
    diagram between them.
"""

# Imports
import union_dataset_function
import sqlite3
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import group_factor


# Function
def get_exons_list(cnx, sf_list, regulation):
    """
    Return every non-redundant exons regulated by  at least one factor in ``sf_list`` (with the regulation \
    ``regulation``)

    :param cnx: (sqlite3 connexion object) allow connexion to sed database
    :param sf_list:  (list of string) the list of splicing factor studied
    :param regulation: (list of string) up or down or up + down
    :return: (list of list of int) list of exons shownig the regulation ``regulation`` at least for a factor \
    in ``sf_list``
    """
    exon_list = []
    for sf_name in sf_list:
        exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    exon_list = union_dataset_function.washing_events_all(exon_list)
    return exon_list


def venn_diagram_creator(list1, name1, list2, name2, output):
    """
    Create a venn diagram of the `list1` and `list2` lists of values.

    :param list1: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name1: (string) the name of the list of exons named `list1`
    :param list2: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name2: (string) the name of the list of exons named `list2`
    :param output: (string) path where the result venn diagram will be created
    :return:
    """
    flat_list1 = ["%s_%s" % (exon[0], exon[1]) for exon in list1]
    flat_list2 = ["%s_%s" % (exon[0], exon[1]) for exon in list2]
    # plt.figure(figsize=(48. / 2.54, 27 / 2.54))
    venn2([set(flat_list1), set(flat_list2)], set_labels=(name1, name2))
    plt.savefig("%sVenn_%s_vs_%s.pdf" % (output, name1, name2))
    plt.clf()
    plt.cla()
    plt.close()


def intersection_supresser_and_file_writer(list1, name1, list2, name2, output):
    """
    Write the list of exons `list1` and `list2` into distinct file and without the common exon \
    in `list1` and `list2`

    :param list1: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name1: (string) the name of the list of exons named `list1`
    :param list2: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name2: (string) the name of the list of exons named `list2`
    :param output: (string) path where the result venn diagram will be created
    """
    new_list1 = [exon for exon in list1 if exon not in list2]
    new_list2 = [exon for exon in list2 if exon not in list1]
    file_writer(new_list1, name1, output)
    file_writer(new_list2, name2, output)


def file_writer(list_exon, name_list, output):
    """
    Write a file containing exons (on per line) identified by their gene_id and their position within the hosting gene.

    :param list_exon: (list of 2 int) list of exons identified by their gene_id and exon_position
    within the hosting gene
    :param name_list: (string) the name of the list of exons named `list_exon`
    :param output: (string) path where the file will be created
    :return:
    """
    with open("%s%s_exons" % (output, name_list), "w") as my_file:
        for exon in list_exon:
            my_file.write("%s\t%s\n" % (exon[0], exon[1]))


def main():

    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/")
    div_group = {"AT_rich": group_factor.at_rich_down, "GC_rich": group_factor.gc_rich_down,
                 "AT_rich_U2": list(group_factor.at_rich_down) + list(group_factor.u2_factors),
                 "GC_rich_U1": list(group_factor.gc_rich_down) + list(group_factor.u1_factors),
                 "U1": group_factor.u1_factors, "U2": group_factor.u2_factors}
    dic_exon = {}
    for name_group in div_group.keys():
        print("Getting all exon regulated by %s factor" % name_group)
        dic_exon[name_group] = get_exons_list(cnx, div_group[name_group], "down")
    for key in dic_exon.keys():
        file_writer(dic_exon[key], "%s_with_intersection" % key, output)
    for couple in [["AT_rich", "GC_rich"], ["AT_rich_U2", "GC_rich_U1"], ["U1", "U2"]]:
        venn_diagram_creator(dic_exon[couple[0]], couple[0], dic_exon[couple[1]], couple[1], output)
        intersection_supresser_and_file_writer(dic_exon[couple[0]], couple[0], dic_exon[couple[1]], couple[1], output)
    cnx.close()


if __name__ == "__main__":
    main()
