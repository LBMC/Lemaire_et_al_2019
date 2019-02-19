#!/usr/bin/python3

# -*- coding: utf-8 -*-
import os
import exon_control_handler
import numpy as np
import figure_producer
import sys
import control_bp_ppt
import control_mfe


def get_control_exon_size_information(cnx, exon_type):
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
        query = """SELECT upstream_intron_size, exon_size, downstream_intron_size,
                   iupac_upstream_intron, iupac_downstream_intron
                   FROM sed
                   WHERE exon_type LIKE '%{}%'""".format(exon_type)
    else:
        query = """SELECT upstream_intron_size, exon_size, downstream_intron_size,
                   iupac_upstream_intron, iupac_downstream_intron
                   FROM sed
                """
    cursor.execute(query)
    names = [description[0] for description in cursor.description]
    result = cursor.fetchall()
    # turn tuple into list
    nresult = []
    for exon in result:
        nresult.append(list(exon))
    return nresult


def get_median_iupac_introns(upstream, downstream):
    """
    Get the median value for every iupac nucleotide frequencies in upstream and downstream intron of an exon.

    :param upstream: (string) frequency of each nucleotides/dnt in full upstream intron separated by a ';',
    the last value corresponds to the size of the upstream sequence
    :param downstream: (string) frequency of each nucleotides/dnt in full downstream intron separated by a ';',
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
               "min_flanking_intron_size": [], "iupac_mean_intron": {nt: [] for nt in figure_producer.nt_dic},
               }
    count = 0
    exon_info_len = len(exon_info)
    for exon in exon_info:
        intron = exon[-2:]
        exon = np.array(exon[:-2], dtype=float)
        mean_intron_iupac = get_median_iupac_introns(intron[0], intron[1])
        for key in new_dic["iupac_mean_intron"].keys():
            new_dic["iupac_mean_intron"][key].append(mean_intron_iupac[figure_producer.nt_dic[key]])

        new_dic['median_flanking_intron_size'].append(np.nanmedian([exon[0], exon[2]]))
        new_dic["min_flanking_intron_size"].append(np.nanmin([exon[0], exon[2]]))
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


def get_summary_dictionaries(exons_dictionary):
    """
    Summary a dictionary that contains the data about individual exons.

    :param exons_dictionary: (dictionary of list of values) contains every information about exon in exon_info
    :return: (dictionary) summarized mean sd and median for every exon within ``exons_dictionary``
    """
    median_dic = {}
    for key in exons_dictionary:
        if "iupac" in key:
            median_dic[key] = {}
            if "iupac" in key:
                for nt in exons_dictionary[key].keys():
                    cur_list = np.array(exons_dictionary[key][nt], dtype=float)
                    median_dic[key][nt] = np.median(cur_list[~np.isnan(cur_list)])
        else:
            cur_list = np.array(exons_dictionary[key], dtype=float)
            median_dic[key] = np.median(cur_list[~np.isnan(cur_list)])
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


def control_handler(cnx, exon_type, size=None):
    my_path = os.path.dirname(os.path.realpath(__file__))
    control_folder = my_path + "/control"
    control_file = control_folder + "/control.py"
    control_full = control_folder + "/control_full.pkl"
    ctrl_list, tmp = exon_control_handler.get_control_information(exon_type, control_file, control_full)
    if ctrl_list is None:
        print("Control dictionary was not found !")
        print("Creating control information")
        names, exon_tuple = exon_control_handler.get_control_exon_information(cnx, exon_type)
        # getting the new columns
        exon_tuple = exon_control_handler.remove_redundant_gene_information(exon_tuple)
        tmp = exon_control_handler.create_a_temporary_dictionary(names, exon_tuple)
        ctrl_list = exon_control_handler.get_summary_dictionaries(names, tmp)
        exon_control_handler.write_control_file(exon_type, control_file, str(ctrl_list))
        exon_control_handler.write_pickle(control_full, tmp)
    if "rel_exon_intron_up" not in ctrl_list.keys():
        print("relative exon_intron size where not found.")
        print("getting relative exon_intron size")
        exon_tuple = get_control_exon_size_information(cnx, exon_type)
        print("Relative size calculation...")
        tmp_dic = tmp_dic_creator(exon_tuple)
        tmp = dict(tmp, **tmp_dic)
        print("summarizing")
        sum_dic = get_summary_dictionaries(tmp_dic)
        ctrl_list = dict(ctrl_list, **sum_dic)
        print("writting...")
        write_adapted_dic(control_file, exon_type, str(ctrl_list))
        exon_control_handler.write_pickle(control_full, tmp)
    if size is not None and "nb_good_bp_%s" % size not in ctrl_list.keys():
        file_name = "%s_%s_bp_ppt_score.py" % (exon_type, size)
        ctrl_file = control_folder + "/" + file_name
        if os.path.isfile(ctrl_file):
            print("Loading %s file to add nb_ggod_pb in control dic !" % ctrl_file)
            sys.path.insert(0, control_folder)
            mod = __import__(file_name.replace(".py", ""))
            val = np.median(mod.nb_good_bp)
            ctrl_list["nb_good_bp_%s" % size] = val
            tmp["nb_good_bp_%s" % size] = mod.nb_good_bp
            val = np.median(mod.ag_count)
            ctrl_list["ag_count"] = val
            tmp["ag_count"] = mod.ag_count
            val = np.median(mod.hbound)
            ctrl_list["hbound"] = val
            tmp["hbound"] = mod.hbound
            write_adapted_dic(control_file, exon_type, str(ctrl_list))
            exon_control_handler.write_pickle(control_full, tmp)
        else:
            print("Creating control file, this operation will take some time")
            control_bp_ppt.control_dictionaries_creator(exon_type, size)
    if "mfe_3ss" not in ctrl_list.keys():
        file_name = "%s_mfe.py" % exon_type
        ctrl_file = control_folder + "/" + file_name
        if os.path.isfile(ctrl_file):
            print("Loading %s file to add mfe_3ss and mfe_5ss in control dic !" % ctrl_file)
            sys.path.insert(0, control_folder)
            module = __import__(file_name.replace(".py", ""))
            ctrl_list["mfe_3ss"] = np.median(module.mfe_3ss)
            tmp["mfe_3ss"] = module.mfe_3ss
            ctrl_list["mfe_5ss"] = np.median(module.mfe_5ss)
            tmp["mfe_5ss"] = module.mfe_5ss
            write_adapted_dic(control_file, exon_type, str(ctrl_list))
            exon_control_handler.write_pickle(control_full, tmp)
        else:
            print("Creating control file, for mfe.")
            control_mfe.control_dictionaries_creator()

    return ctrl_list, tmp


if __name__ == "__main__":
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx_sed = figure_producer.connexion(seddb)
    control_handler(cnx_sed, "CCE", None)
