#!/usr/bin/python3.5

# coding: utf-8

import os
import sys
import numpy as np


# function
def get_control_information(exon_type, control_file):
    """
    Test whether we have all the information we need for ``exon_type`` exons

    :param exon_type: (string) the type of control exon we want to use
    :param control_file: (string) the file were we want to write the control data results
    :return: None if the control file does'nt exist, or a dictionary that contains \
    all the information about the CCE exon
    """
    control_folder = os.path.dirname(control_file)
    if not os.path.isdir(control_folder):
        os.mkdir(control_folder)
        return None
    if not os.path.isfile(control_file):
        return None
    sys.path.insert(0, control_folder)
    mod = __import__("control")
    try:
        dic = eval("mod.%s_dic" % exon_type)
    except AttributeError:
        return None
    return dic


def get_control_exon_information(cnx, exon_type):
    """
    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return:
        * result: (list of tuple) every information about control exons
        * names: (list of string) the name of every column in sed table
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT *
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT *
                   FROM sed
                """
    cursor.execute(query)
    names = [description[0] for description in cursor.description]
    result = cursor.fetchall()
    return names, result


def get_iupac_nt(iupac_string, nt):
    """
    Get the frequency of a ``nt`` in a iupac string

    :param nt: (string) the nt we want to find
    :param iupac_string: (string) contains iupac frequency for an exon in that order :  \
    ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"].
    :return: (float) the ``nt`` frequency in a iupac string
    """
    iupac_list = ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]
    iupac_string = iupac_string.replace("\n", "").split(";")
    for i in range(len(iupac_list)):
        if iupac_list[i] == nt:
            return float(iupac_string[i])


def create_a_temporary_dictionary(names, exon_info):
    """
    :param names: (list of string) the name of every column in sed table
    :param exon_info: (list of tuple) each tuple corresponds to an exon \
    (or a row in sed table)
    :return: (dictionary of list of values) contains every information about exon in exon_info
    """
    iupac_list = ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]
    new_dic = {}
    # initialization of the new dictionary
    for key in names:
        if "iupac" not in key:
            new_dic[key] = []
        else:
            new_dic[key] = {nt: [] for nt in iupac_list}

    for i in range(len(exon_info)):
        for j in range(len(exon_info[i])):
            if "iupac" not in names[j]:
                new_dic[names[j]].append(exon_info[i][j])
            else:
                if exon_info[i][j] is None:
                    for nt in iupac_list:
                        new_dic[names[j]][nt].append(None)
                else:
                    for nt in iupac_list:
                        new_dic[names[j]][nt].append(get_iupac_nt(exon_info[i][j], nt))
    return new_dic


def get_summary_dictionaries(names, exons_dictionary):
    """
    Summary a dictionary that contains the data about individual exons.
    :param names: (list of string) the name of every column in sed table
    :param exons_dictionary: (dictionary of list of values) contains every information about exon in exon_info
    :return: (dictionary) summarized mean sd and median for every exon within ``exons_dictionary``
    """
    median_dic = {}
    # initialization
    for key in names[4::]:
        if "iupac" not in key:
            if exons_dictionary[key]:
                cur_list = np.array(exons_dictionary[key], dtype=float)
                median_dic[key] = np.median(cur_list[~np.isnan(cur_list)])
            else:
                median_dic[key] = [float("nan")]
        else:
            median_dic[key] = {}
            for nt in ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]:
                if exons_dictionary[key][nt]:
                    cur_list = np.array(exons_dictionary[key][nt], dtype=float)
                    median_dic[key][nt] = np.median(cur_list[~np.isnan(cur_list)])
                else:
                    median_dic[key][nt] = [float("nan")]
    return median_dic


def write_control_file(exon_type, control_file, str2write):
    """
    :param exon_type: (string) the type of control exon we want to use
    :param control_file: (string) the file were we want to write the control data results
    :param str2write:  (string) the result we want to write in the file
    """
    with open(control_file, "a") as ctrl_file:
        ctrl_file.write("%s_dic = %s\n" % (exon_type, str2write))


def control_handler(cnx, exon_type):
    my_path = os.path.dirname(os.path.realpath(__file__))
    control_folder = my_path + "/control"
    control_file = control_folder + "/control.py"
    ctrl_dic = get_control_information(exon_type, control_file)
    if ctrl_dic is None:
        print("Control dictionary was not found !")
        print("Creating control information")
        names, exon_tuple = get_control_exon_information(cnx, exon_type)
        tmp = create_a_temporary_dictionary(names, exon_tuple)
        ctrl_dic = get_summary_dictionaries(names, tmp)
        write_control_file(exon_type, control_file, str(ctrl_dic))
    return ctrl_dic
