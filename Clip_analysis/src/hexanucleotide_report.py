#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    The goal of this script is from a result folder produced by launcher.py script, to get for every splicing factor\
    the N most frequent hexanucleotides found in their clip experiment and build weblogo of those hexanucleotide.
"""

import subprocess
import os
import pandas as pd
import figure_maker
import argparse


def frequencie_file_finder(folder):
    """
    Find every frequency files found in the folder ``folder``.

    :param folder: (string) folder containing the hexanucleotide frequencies file.
    :return: (dictionary of list of string) each key of the dictionary is link with file having the same basename.
    """
    dic_files = {}
    list_file = subprocess.check_output("find %s -name \"frequencies*\" -type f" % folder,
                                        shell=True, stderr=subprocess.STDOUT).decode("ascii")
    list_file = list_file.split("\n")[:-1]
    for my_file in list_file:
        base = os.path.basename(my_file).split(".")[0]
        if base not in dic_files.keys():
            dic_files[base] = [my_file]
        else:
            dic_files[base].append(my_file)
    return dic_files


def get_splicing_factor_studied(list_files):
    """
    Get the name of the factors studied from every file in ``list_files``.

    :param list_files: (lits of string) list of hexanucleotide frequency file
    :return: (list of string) the list of factors studied
    """
    sf_list = []
    for myfile in list_files:
        sf_name = myfile.split("/")[-3].split("_")[0].upper()
        if sf_name not in sf_list:
            sf_list.append(sf_name)
    return sf_list


def get_sum_frequencies(list_files):
    """
    Get the average frequencies of every hexanucleotide of ``list_files``.

    :param list_files: (list of string) list of files
    :return: (pandas DataFrame) the average frequencies of every hexanucleotide found in ``list_files``
    """
    df_sum = pd.read_csv(list_files[0], sep="\t", index_col=0)
    if len(list_files) == 1:
        return df_sum.sort_values("frequencies", ascending=False)
    for i in range(1, len(list_files)):
        df = pd.read_csv(list_files[i], sep="\t", index_col=0)
        df_sum += df
    df_sum["frequencies"] = df_sum["frequencies"] / len(list_files)
    return df_sum.sort_values("frequencies", ascending=False)


def result_writer(list_files, output, name_template, nb_seq):
    """
    Write the list of the ``nb_seq`` most frequent hexanucleotides and make a weblogo of them.

    :param list_files: (list of string) the list of hexanucleotide frequency files in the clip of a factor \
    in a given template.
    :param output:  (string) path where the result will be created
    :param name_template: (string) the name of the template used
    :param nb_seq:  (int) the number of hexanucleotides selected
    """
    filename = "%s/hexanucleotide_frequencies_%s.txt" % (output, name_template)
    print(filename)
    df = get_sum_frequencies(list_files)
    df = df.head(nb_seq)
    df.to_csv(filename, sep="\t")
    hexant = list(df.index)
    print(hexant)
    sf_name = output.split("/")[-1]
    if "intron" in name_template:
        web_name = "Clip_%s_intron" % sf_name
    elif "exon" in name_template:
        web_name = "Clip_%s_exon" % sf_name
    else:
        web_name = "Clip_%s_all" % sf_name
    figure_maker.web_logo_creator(hexant, web_name, output + "/")


def main(folder, output, nb_seq):
    """
    Create the weblogo of the most enriched hexanucleotides.

    :param folder: (string) folder containing the result frequency files
    :param output: (string) folder were the result will be created
    :param nb_seq: (int) the number of hexanucleotide to consider
    """
    dic_files = frequencie_file_finder(folder)
    for key in dic_files.keys():
        list_files = dic_files[key]
        sf_list = get_splicing_factor_studied(list_files)
        for sfname in sf_list:
            new_list_file = []
            for my_file in list_files:
                if sfname in my_file:
                    new_list_file.append(my_file)
            if len(new_list_file) > 0:
                my_folder = output + "/" + sfname
                if not os.path.isdir(my_folder):
                    os.mkdir(my_folder)
                print(my_folder)
                result_writer(new_list_file, my_folder, key, nb_seq)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    From a folder containing the hexanucleotide frequency files produced by the script launcher.py.
    Get for each splicing factor the hexanucleotide having the highest frequencies in clip data, and create a weblogo
    of them.
    """)
    # Arguments for the parser

    req_arg = parser.add_argument_group("Required arguments")

    req_arg.add_argument('--folder', dest='folder', help="a folder containing the frequency files",
                         required=True)
    req_arg.add_argument('--output', dest='output', help="a folder where the result will be created",
                         required=True)
    parser.add_argument('--nb_seq', dest='nb_seq',
                        help="The number of the most frequence hexanucleotide to consider", default=10)

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        parser.error("The output directory doesn't exists !")
    if not os.path.isdir(args.folder):
        parser.error("The input folder doe not exists !")
    try:
        args.nb_seq = int(args.nb_seq)
    except ValueError:
        parser.error("Wrong value for the argument --nb_seq, it must be an integer")

    if args.output[-1] == "/":
        args.output = args.output[:-1]

    main(args.folder, args.output, args.nb_seq)


if __name__ == "__main__":
    launcher()
