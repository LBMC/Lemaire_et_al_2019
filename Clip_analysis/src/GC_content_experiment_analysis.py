#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is to create 3 table that show the gc content of every clip project analyzed
"""

# IMPORTS
import pandas as pd
import subprocess
import numpy as np
import os
import sys



# FUNCTIONS
def get_template_name(folder):
    """
    Gave the name of every frequencies files located in a subdirectory frequencies within a folder containing every \
    clip seq result.

    :param folder: (string) folder that contains every clip seq results
    :return: (list of strings) the name of every result folder
    """
    cmd = "ls %s/*/frequencies/*.py" % folder
    list_files = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode("ascii")
    list_files = list_files.split("\n")[:-1]
    new_files = []
    for my_file in list_files:
        new_files.append(my_file.split("/")[-1])
    return np.unique(new_files)


def get_interest_list_files(folder, template):
    """
    Get the files named ``template`` in the folder ``folder``.

    :param folder:  (string) folder that contains every clip seq results
    :param template: (string) the names of the files of interest
    :return: (list of string) list of every files named ``template`` located within ``folder``
    """
    cmd = "find %s -name %s -type f" % (folder, template)
    list_files = subprocess.check_output(cmd.split(" "), stderr=subprocess.STDOUT).decode("ascii")
    return list_files.split("\n")[:-1]


def load_gc_content_values(file_name):
    """
    Get the gc content value located in the file ``file_name``.

    :param file_name: (string) python files containing a dictionary of nucleotide frequencies
    :return: (list of one string + one float) the name of the clip project and it's associated gc frequencies
    """
    project = file_name.split("/")[-3]
    with open(file_name) as outfile:
        for line in outfile:
            line = line.split(":")
            for i in range(len(line)):
                if "S" in line[i]:
                    freq = line[i+1].split(",")[0]
                    freq = round(float(freq), 2)
                    return [project, str(freq).replace(".", ",")]


def load_gc_content_values_from_list(res_dic, list_files):
    """
    Get the gc content value located in every file of ``list_files``.

    :param res_dic: (dictionary of list of float) link each project (key) to the frequencies of gc \
    in each one of the template (value)
    :param list_files: (list of string) list of interest files from which we want to recover the gc content.
    :return: (dictionary of list of float) res dic updated
    """
    for my_file in list_files:
        gc_info = load_gc_content_values(my_file)
        if gc_info[0] in res_dic.keys():
            res_dic[gc_info[0]].append(gc_info[1])
        else:
            res_dic[gc_info[0]] = [gc_info[1]]
    return res_dic


def get_result_dic(folder, list_template):
    """
    Get the gc content for every template of every clip seq experiment.

    :param folder:  (string) folder that contains every clip seq results
    :param list_template: (list of string) list of template files
    :return:
        res_dic : (dictionary of list of float) (dictionary of list of float) link each project (key) to the frequencies of gc \
        in each one of the template (value)
        index : (list of strings) the name of the template
    """
    res_dic = {}
    index = []
    for template in list_template:
        print("########## TEMPLATE ###################")
        index.append(template.replace("frequencies_", ""))
        list_files = get_interest_list_files(folder, template)
        res_dic = load_gc_content_values_from_list(res_dic, list_files)
    return res_dic, index


def create_result_table(res_dic, index):
    """
    Create a pandas dataframe from ``res_dic``and ``index``
    :param res_dic: (dictionary of list of float) (dictionary of list of float) link each project (key) to the frequencies of gc \
        in each one of the template (value)
    :param index: (list of strings) the name of the template
    :return: (pandas DataFrame) table of GC content
    """
    df = pd.DataFrame(res_dic, index=index)
    return df.transpose()


def main(folder):
    """
    Create the GC frequencies table of every clip project, for every template used for the clip analysis

    :param folder: (string) folder where the gc content of every clip seq project are located
    """
    output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/gc_frequencies_report/")
    if not os.path.isdir(output):
        os.mkdir(output)
    list_template = get_template_name(folder)
    res_dic, index = get_result_dic(folder, list_template)
    df = create_result_table(res_dic, index)
    table_name = "%sgc_content_analysis.txt" % output
    df.to_csv(table_name, sep="\t")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main("/media/nicolas/DD_1/Splicing_Lore_project/Clip_analysis/result/GC_content_analysis/")