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
from matplotlib_venn import venn2, venn3
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


def venn_diagram_creator(list1, name1, color1, list2, name2, color2, output):
    """
    Create a venn diagram of the `list1` and `list2` lists of values.

    :param list1: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name1: (string) the name of the list of exons named `list1`
    :param color1: (string) the color of list1
    :param list2: (list of 2 int) list of exons identified by their gene_id and exon_position within the hosting gene
    :param name2: (string) the name of the list of exons named `list2`
    :param color2: (string) the color of list 2
    :param output: (string) path where the result venn diagram will be created
    :return:
    """
    flat_list1 = ["%s_%s" % (exon[0], exon[1]) for exon in list1]
    flat_list2 = ["%s_%s" % (exon[0], exon[1]) for exon in list2]
    # plt.figure(figsize=(48. / 2.54, 27 / 2.54))
    c = venn2([set(flat_list1), set(flat_list2)], set_labels=(name1, name2))
    c.get_patch_by_id('10').set_color(color1)
    c.get_patch_by_id('01').set_color(color2)
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

    :param list_exon: (list of 2 int) list of exons identified by their gene_id and exon_position \
    within the hosting gene
    :param name_list: (string) the name of the list of exons named `list_exon`
    :param output: (string) path where the file will be created
    """
    with open("%s%s_exons" % (output, name_list), "w") as my_file:
        for exon in list_exon:
            my_file.write("%s\t%s\n" % (exon[0], exon[1]))


def venn3_diagram_creator(list_exons, list_names, list_color, output):
    """
    Create a venn diagram of the `list1` and `list2` lists of values.

    :param list_exons: (list of list of 2 int) list of exons
    :param list_names: (list of string) the names of the list of exons
    :param output: (string) path where the result venn diagram will be created
    :return:
    """
    def flaten_list(x): return ["%s_%s" % (a[0], a[1]) for a in x]
    new_lists = [set(flaten_list(my_list)) for my_list in list_exons]
    c = venn3(new_lists, set_labels=list_names)
    c.get_patch_by_id('100').set_color(list_color[0])
    c.get_patch_by_id('010').set_color(list_color[1])
    c.get_patch_by_id('001').set_color(list_color[2])
    plt.savefig("%sVenn_%s.pdf" % (output, "_".join(list_names)))
    plt.clf()
    plt.cla()
    plt.close()


def main():
    color_dic2 = {"GC-exons": "#5555FF", "AT-exons": "#00aa00",
                 "U1-exons": "#7777F0", "U2-exons": "#33EE33",
                 "AT-exons_U2": "#aaFFaa", "GC_rich_U1": "#aaaaFF"}
    couples = [["AT-exons", "GC-exons"], ["U1-exons", "U2-exons"], ["GC-exons", "U1-exons"], ["AT-exons", "U1-exons"], ["GC-exons", "U2-exons"], ["AT-exons", "U2-exons"]]
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/")
    div_group = {"AT-exons-all": group_factor.at_rich_down, "GC-exons-all": group_factor.gc_rich_down,
                 "AT-exons_U2": list(group_factor.at_rich_down) + list(group_factor.u2_factors),
                 "GC-exons_U1": list(group_factor.gc_rich_down) + list(group_factor.u1_factors),
                 "U1-exons": ["SNRNP70", "SNRPC", "DDX5_DDX17"],
                 "U2-exons": ["SF1", "SF3A3", "SF3B4", "U2AF2"]}
    dic_exon = {}
    for name_group in div_group.keys():
        exon_list = get_exons_list(cnx, div_group[name_group], "down")
        print("Exon regulated by %s : %s" % (name_group, len(exon_list)))
        dic_exon[name_group] = exon_list
    # for key in dic_exon.keys():
    #     file_writer(dic_exon[key], "%s_with_intersection" % key, output)
    dic_exon["GC-exons"] = [exon for exon in dic_exon["GC-exons-all"]
                            if exon not in dic_exon["AT-exons-all"]]
    dic_exon["AT-exons"] = [exon for exon in dic_exon["AT-exons-all"]
                            if exon not in dic_exon["GC-exons-all"]]
    print("GC-exons : %s exons" % len(dic_exon["GC-exons"]))
    print("AT-exons : %s exons" % len(dic_exon["AT-exons"]))
    for couple in couples:
        venn_diagram_creator(dic_exon[couple[0]], couple[0], color_dic2[couple[0]],
                             dic_exon[couple[1]], couple[1], color_dic2[couple[1]], output)
        intersection_supresser_and_file_writer(dic_exon[couple[0]], couple[0],
                                               dic_exon[couple[1]], couple[1],
                                               output)
    list_name1 = ["U1-exons", "AT-exons", "GC-exons"]
    list_color = [color_dic2[mylist] for mylist in list_name1]
    list_exon = [dic_exon[mylist] for mylist in list_name1]
    venn3_diagram_creator(list_exon, list_name1, list_color, output)
    list_name2 = ["U2-exons", "AT-exons", "GC-exons"]
    list_color = [color_dic2[mylist] for mylist in list_name2]
    list_exon = [dic_exon[mylist] for mylist in list_name2]
    venn3_diagram_creator(list_exon, list_name2, list_color, output)
    cnx.close()


if __name__ == "__main__":
    main()
