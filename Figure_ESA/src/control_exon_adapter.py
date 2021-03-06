#!/usr/bin/python3

# -*- coding: utf-8 -*-
import os
import exon_control_handler
import numpy as np
import functions
import union_dataset_function
import sys


def get_control_exon_size_information(cnx, exon_type, exon2remove):
    """
    Get the surrounding intron and the exon size information.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos, upstream_intron_size, exon_size, downstream_intron_size,
                   iupac_upstream_intron, iupac_downstream_intron
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT gene_id, exon_pos, upstream_intron_size, exon_size, downstream_intron_size,
                   iupac_upstream_intron, iupac_downstream_intron
                   FROM sed
                """
    cursor.execute(query)
    names = [description[0] for description in cursor.description]
    result = cursor.fetchall()
    # turn tuple into list
    print("in function get_control_exon_size_information, CCE exon before removing bad ones : %s" % len(result))
    nresult = [list(exon)[2:] for exon in result if [exon[0], exon[1]] not in exon2remove]
    print("CCE exon after removing those regulated by splicing factors : %s" % len(nresult))
    return nresult


def get_median_iupac_introns(upstream, downstream):
    """
    Get the median value for every iupac nucleotide frequencies in upstream and downstream intron of an exon.

    :param upstream: (string) frequency of each nucleotides/dnt in full upstream intron separated by a ';' \
    the last value corresponds to the size of the upstream sequence
    :param downstream: (string) frequency of each nucleotides/dnt in full downstream intron separated by a ';' \
    the last value corresponds to the size of the downstream sequence
    :return: (list of float) the frequencies of nt and dnt in the introns
    """
    if upstream is not None and downstream is not None:
        upstream = list(map(float, upstream.split(";")))
        downstream = list(map(float, downstream.split(";")))
        result = []
        for i in range(len(upstream)):
            result.append(np.mean([upstream[i], downstream[i]]))
        return result
    else:
        return [None] * 10


def tmp_dic_creator(exon_info):
    """
    :param exon_info: (list of tuple) each tuple corresponds to an exon \
    (or a row in sed table)
    :return: (dictionary of list of values) contains relative_size info information about exon in exon_info
    """
    new_dic = {"rel_exon_intron_up": [], "rel_exon_intron_down": [], "rel_exon_introns": [],
               'median_flanking_intron_size': [],
               "min_flanking_intron_size": [], "iupac_mean_intron": {nt: [] for nt in functions.nt_dic},
               "iupac_introns": {nt: [] for nt in functions.nt_dic},
               "introns_size": []}
    count = 0
    exon_info_len = len(exon_info)
    for exon in exon_info:
        intron = exon[-2:]
        exon = np.array(exon[:-2], dtype=float)
        mean_intron_iupac = get_median_iupac_introns(intron[0], intron[1])
        if intron[0] is not None:
            upstream = list(map(float, intron[0].split(";")))
        else:
            upstream = [None] * 10
        if intron[1] is not None:
            downstream = list(map(float, intron[1].split(";")))
        else:
            downstream = [None] * 10
        for key in new_dic["iupac_mean_intron"].keys():
            new_dic["iupac_mean_intron"][key].append(mean_intron_iupac[functions.nt_dic[key]])
            new_dic["iupac_introns"][key].append(upstream[functions.nt_dic[key]])
            new_dic["iupac_introns"][key].append(downstream[functions.nt_dic[key]])

        new_dic['median_flanking_intron_size'].append(np.nanmedian([exon[0], exon[2]]))
        new_dic["min_flanking_intron_size"].append(np.nanmin([exon[0], exon[2]]))
        new_dic["introns_size"].append(exon[0])
        new_dic["introns_size"].append(exon[2])
        if exon[0] != 0:
            new_dic["rel_exon_intron_up"].append(float(exon[1]) / float(exon[0]))
        else:
            new_dic["rel_exon_intron_up"].append(None)
        if exon[2] != 0:
            new_dic["rel_exon_intron_down"].append(float(exon[1]) / float(exon[2]))
        else:
            new_dic["rel_exon_intron_down"].append(None)
        if exon[0] != 0 or exon[1] != 0:
            mean = float(exon[0] + exon[2]) / 2
            new_dic["rel_exon_introns"].append(float(exon[1]) / mean)
        else:
            new_dic["rel_exon_introns"].append(None)
        count += 1
        sys.stdout.write("%s / %s \r" % (count, exon_info_len))
        sys.stdout.flush()
    return new_dic


def get_summary_dictionaries(exons_dictionary, summary):
    """
    Summary a dictionary that contains the data about individual exons.

    :param exons_dictionary: (dictionary of list of values) contains every information about exon in exon_info
    :param summary: (string) mean or median
    :return: (dictionary) summarized mean sd and median for every exon within ``exons_dictionary``
    """
    median_dic = {}
    for key in exons_dictionary:
        if "iupac" in key:
            median_dic[key] = {}
            if "iupac" in key:
                for nt in exons_dictionary[key].keys():
                    cur_list = np.array(exons_dictionary[key][nt], dtype=float)
                    median_dic[key][nt] = eval("np.%s(cur_list[~np.isnan(cur_list)])" % summary)
        else:
            cur_list = np.array(exons_dictionary[key], dtype=float)
            median_dic[key] = eval("np.%s(cur_list[~np.isnan(cur_list)])" % summary)
    return median_dic


def write_adapted_dic(control_file, exon_type, str2write):
    """
    Write the dictionary ``str2write`` in the ``control_file``.

    :param control_file: (string) the control file
    :param exon_type: (string) the control exon type
    :param str2write: ( the string to write) the doc 2 write
    """
    list_of_lines = []
    with open(control_file, 'r') as outfile:
        line = outfile.readline()
        while line:
            list_of_lines.append(line)
            line = outfile.readline()
    with open(control_file, 'w') as outfile:
        for line in list_of_lines:
            if "%s_dic" % exon_type in line:
                outfile.write("%s_dic = %s\n" % (exon_type, str2write))
            else:
                outfile.write("%s\n" % line)


def control_handler(cnx, exon_type, summary, regulation="down"):
    my_path = os.path.dirname(os.path.realpath(__file__))
    control_folder = my_path + "/control"
    control_file = "%s/control_%s.py" % (control_folder, summary)
    control_full = control_folder + "/control_full.pkl"
    ctrl_list, tmp = exon_control_handler.get_control_information(exon_type, control_file, control_full)
    exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx, regulation)
    if ctrl_list is None:
        print("Control dictionary was not found !")
        print("Creating control information")
        names, exon_tuple = exon_control_handler.get_control_exon_information(cnx, exon_type, exon2remove, regulation)
        # getting the new columns
        exon_tuple = exon_control_handler.remove_redundant_gene_information(exon_tuple)
        tmp = exon_control_handler.create_a_temporary_dictionary(names, exon_tuple)
        ctrl_list = exon_control_handler.get_summary_dictionaries(names, tmp, summary)
        exon_control_handler.write_control_file(exon_type, control_file, str(ctrl_list))
        exon_control_handler.write_pickle(control_full, tmp)
    if "rel_exon_intron_up" not in ctrl_list.keys():
        print("relative exon_intron size where not found.")
        print("getting relative exon_intron size")
        exon_tuple = get_control_exon_size_information(cnx, exon_type, exon2remove)
        print("Relative size calculation...")
        tmp_dic = tmp_dic_creator(exon_tuple)
        tmp = dict(tmp, **tmp_dic)
        print("summarizing")
        sum_dic = get_summary_dictionaries(tmp_dic, summary)
        ctrl_list = dict(ctrl_list, **sum_dic)
        print("writting...")
        write_adapted_dic(control_file, exon_type, str(ctrl_list))
        exon_control_handler.write_pickle(control_full, tmp)
    return ctrl_list, tmp


if __name__ == "__main__":
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx_sed = functions.connexion(seddb)
    control_handler(cnx_sed, "CCE", "median")
