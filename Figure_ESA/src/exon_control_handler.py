#!/usr/bin/python3.5

# coding: utf-8

import os
import sys
import numpy as np
import pickle
import union_dataset_function


# Function
def get_control_information(exon_type, control_file, control_full):
    """
    Test whether we have all the information we need for ``exon_type`` exons

    :param exon_type: (string) the type of control exon we want to use
    :param control_file: (string) the file were we want to write the control data results
    :param control_full: (string) the file were we want to write the full control data results
    :return: None if the control file does'nt exist, or a dictionary that contains \
    all the information about the CCE exon
    """
    control_folder = os.path.dirname(control_file)
    if not os.path.isdir(control_folder):
        os.mkdir(control_folder)
        return None, None
    if not os.path.isfile(control_file):
        return None, None
    sys.path.insert(0, control_folder)
    mod = __import__(os.path.basename(control_file).replace(".py", ""))
    try:
        dic = eval("mod.%s_dic" % exon_type)
    except AttributeError:
        return None, None
    tmp = load_pickle(control_full)
    return dic, tmp


def get_control_exon_information(cnx, exon_type, exon_2_remove, regulation):
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
    # turn tuple into list
    print("CCE exons : %s" % len(result))
    nresult = [list(exon) for exon in result if [exon[1], exon[2]] not in exon_2_remove]
    print("CCE exons not %s-regulated by splicing factors : %s" % (regulation, len(nresult)))
    return names, nresult


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


def get_dnt(dnt_string, dnt):
    """
    Get the frequency of a ``dnt`` in a dnt string

    :param dnt: (string) the dnt we want to find
    :param dnt_string: (string) contains dnt frequency for an exon in that order :  \
    ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"].
    :return: (float) the ``nt`` frequency in a iupac string
    """
    dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    dnt_string = dnt_string.replace("\n", "").split(";")
    for i in range(len(dnt_list)):
        if dnt_list[i] == dnt:
            return float(dnt_string[i])


def create_a_temporary_dictionary(names, exon_info):
    """
    :param names: (list of string) the name of every column in sed table
    :param exon_info: (list of tuple) each tuple corresponds to an exon \
    (or a row in sed table)
    :return: (dictionary of list of values) contains every information about exon in exon_info
    """
    iupac_list = ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]
    dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    new_dic = {}
    # initialization of the new dictionary
    for key in names:
        if "iupac" in key:
            new_dic[key] = {nt: [] for nt in iupac_list}
        elif "dnt" in key:
            new_dic[key] = {dnt: [] for dnt in dnt_list}
        else:
            new_dic[key] = []

    for i in range(len(exon_info)):
        for j in range(len(names)):
            if "iupac" in names[j]:
                if exon_info[i][j] is None:
                    for nt in iupac_list:
                        new_dic[names[j]][nt].append(None)
                else:
                    for nt in iupac_list:
                        new_dic[names[j]][nt].append(get_iupac_nt(exon_info[i][j], nt))
            elif "dnt" in names[j]:
                if exon_info[i][j] is None:
                    for dnt in dnt_list:
                        new_dic[names[j]][dnt].append(None)
                else:
                    for dnt in dnt_list:
                        new_dic[names[j]][dnt].append(get_dnt(exon_info[i][j], dnt))
            else:
                new_dic[names[j]].append(exon_info[i][j])
    return new_dic


def get_summary_dictionaries(names, exons_dictionary, summary):
    """
    Summary a dictionary that contains the data about individual exons.
    :param names: (list of string) the name of every column in sed table
    :param exons_dictionary: (dictionary of list of values) contains every information about exon in exon_info
    :param summary: (string) indicates the type of summary dictionary we want to create
    :return: (dictionary) summarized mean sd and median for every exon within ``exons_dictionary``
    """
    median_dic = {}
    # initialization
    for key in names[4::]:
        if "iupac" in key or "dnt" in key:
            median_dic[key] = {}
            if "iupac" in key:
                nt_dnt_list = ["A", "C", "G", "T", "S", "W", "R", "Y", "K", "M"]
            else:
                nt_dnt_list = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                               "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
            for nt in nt_dnt_list:
                if exons_dictionary[key][nt]:
                    cur_list = np.array(exons_dictionary[key][nt], dtype=float)
                    median_dic[key][nt] = eval("np.%s(cur_list[~np.isnan(cur_list)])" % summary)
                else:
                    median_dic[key][nt] = [float("nan")]
        else:
            if exons_dictionary[key]:
                cur_list = np.array(exons_dictionary[key], dtype=float)
                median_dic[key] = eval("np.%s(cur_list[~np.isnan(cur_list)])" % summary)
            else:
                median_dic[key] = [float("nan")]
    return median_dic


def write_control_file(exon_type, control_file, str2write, name=True):
    """
    :param exon_type: (string) the type of control exon we want to use
    :param control_file: (string) the file were we want to write the control data results
    :param str2write:  (string) the result we want to write in the file
    :param name: (boolean) True to write the name `exon_type`_dic in the control file, False else.
    """
    with open(control_file, "a") as ctrl_file:
        if name:
            ctrl_file.write("%s_dic = %s\n" % (exon_type, str2write))
        else:
            ctrl_file.write(str2write)


def remove_redundant_gene_information(exon_list):
    """
    Remove redundant gene data.

    If two exons are regulated by the same gene, the information at gene level for those exons \
    will only be used once.

    :param exon_list: (2 list of list (one for up exon the other for down exon)) that contains every exons \
     given in input by the user.
    :return: the input exons list without redundant data at gene level
    """
    gene_dic = {}
    redundant_dic = {}
    for exon in exon_list:
        if exon[1] not in gene_dic:
            gene_dic[exon[1]] = 1
        else:
            gene_dic[exon[1]] += 1
            redundant_dic[exon[1]] = gene_dic[exon[1]]
    for gene in redundant_dic:
        for exon in exon_list:
            if exon[1] == gene:
                exon[4:9] = [None] * 5
                redundant_dic[gene] -= 1
            if redundant_dic[gene] == 1:
                break
    return exon_list


def write_pickle(filename, dictionary):
    """
    Write a pickle dictionary.

    :param filename: (string) the name used to write the pickle dictionary
    :param dictionary: (dictionary) a dictionary
    """
    output = open(filename, 'wb')
    pickle.dump(dictionary, output)
    output.close()


def load_pickle(filename):
    """
    load a pickle dictionary.

    :param filename: (string) the name used to write the pickle dictionary
    :return: (dictionary) the dictionary in filename
    """
    output = open(filename, 'rb')
    my_dic = pickle.load(output)
    output.close()
    return my_dic


def control_handler(cnx, exon_type, regulation, summary="median"):
    my_path = os.path.dirname(os.path.realpath(__file__))
    control_folder = my_path + "/control"
    control_file = "%s/control_%s.py" % (control_folder, summary)
    control_full = control_folder + "/control_full.pkl"
    ctrl_list, tmp = get_control_information(exon_type, control_file, control_full)
    if ctrl_list is None:
        print("Control dictionary was not found !")
        print("Creating control information")
        exon2remove = union_dataset_function.get_exon_regulated_by_sf(cnx, regulation)
        names, exon_tuple = get_control_exon_information(cnx, exon_type, exon2remove, regulation)
        # getting the new columns
        exon_tuple = remove_redundant_gene_information(exon_tuple)
        tmp = create_a_temporary_dictionary(names, exon_tuple)
        ctrl_list = get_summary_dictionaries(names, tmp, summary)
        write_control_file(exon_type, control_file, str(ctrl_list))
        write_pickle(control_full, tmp)
    return ctrl_list, tmp
