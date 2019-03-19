#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import sqlite3
import os
import exon_class_metaexon
import metaexon_figure_creator
import argparse
import win_size


def extract_exon_list(a_file):
    """
    Extract the exon contained in ``a_file``. Each line in ``a_file`` must be like ``gene_id exon_position`` \
    values in this files have to be tab separated.
    :param a_file: (string) path to a file containing and exons
    :return: (list of lists of 2 int) each sublist corresponds to an exons (gene_id and exon position)
    """
    exon_list = []
    with open(a_file, "r") as in_file:
        line = in_file.readline()
        while line:
            line = line.replace("\n", "")
            exon_list.append(list(map(int, line.split("\t"))))
            line = in_file.readline()
    return exon_list


def main(files, name_files, nt_list, name_fig, exon_type, color_list, legend):
    """
    Create meta_exon figures from the list of exons given in ``files``
    :param files: (list of strings) a file or many files, if many file are given, they have to be coma separated.
    :param name_files: (list of strings) name of the files given in ``files`` if many names are given, they have to \
    be coma separated.
    :param nt_list: (list of strings) a nt or a many nt, if many nt are given they must be separated by a coma.
    :param name_fig: (string) the name of the wanted meta-exon
    :param exon_type: (string) the type of control exons we want to use
    :param color_list: (list of strings) the list of color we want to use
    :param legend: (boolean) True if you want to draw an exon + legend false else
    """
    exon_class_metaexon.set_debug(0)
    window_size = win_size.window_size
    output = os.path.realpath(os.path.dirname(__file__)).replace("src/metaexon_figure", "result/metaexon_figure/")
    if not os.path.isdir(output):
        os.mkdir(output)
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src/metaexon_figure", "data/fasterDB_lite.db")
    cnx = sqlite3.connect(fasterdb)
    list_of_vector_5p = []
    list_of_vector_3p = []
    for i in range(len(files)):
        exon_list = extract_exon_list(files[i])
        new_exon_list = [exon_class_metaexon.ExonClass(cnx, None, exon[0], exon[1], window_size) for exon in exon_list]
        final_res_5p, final_res_3p, p5_analyzed, p3_analyzed = \
            exon_class_metaexon.get_metagene_vectors_windowsed(new_exon_list, window_size)
        list_of_vector_5p.append(final_res_5p)
        list_of_vector_3p.append(final_res_3p)

    for nt in nt_list:
        metaexon_figure_creator.make_metagene_graphics(output, list_of_vector_5p, list_of_vector_3p, nt,
                                                       name_fig, name_files, color_list, exon_type, legend)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    This program aims to create a meta exon figures from a list of exon files.
    """)
    # Arguments for the parser

    req_arg = parser.add_argument_group("required arguments")

    req_arg.add_argument('--files', dest='files', help="an input file/or list of files coma separated",
                         required=True)
    req_arg.add_argument('--name_files', dest='name_files',
                         help="The list of name_files corresponding to the files given in files "
                              "(each name must be coma seaprated)",
                         )
    parser.add_argument('--nt', dest='nt', help="nucleotide wanted on the metaexon figures",
                        default="A,C,G,T")
    parser.add_argument('--name_fig', dest='name_fig', help="name of the figure",
                        default="input")
    parser.add_argument('--exon_type', dest='exon_type', help="control exons",
                        default="CCE")
    parser.add_argument('--color', dest='color', help="the color that will be used in the metaexons figures for each"
                        "input given; (coma separated)",
                        default="#0000fe,#00aaab,red,green,purple")
    parser.add_argument('--legend', dest='legend', help="True to display the legend, false else", default=True)

    args = parser.parse_args()

    input_files = args.files.split(",")
    for my_file in input_files:
        if not os.path.isfile(my_file):
            parser.error("The file %s doesn't exist!" % my_file)
    name_files = args.name_files.split(",")
    if len(name_files) != len(input_files):
        parser.error("The number of files and name_files given are different")

    nt_list = args.nt.split(",")
    for nt in nt_list:
        if nt not in ["A", "T", "G", "C", "S", "W", "K", "M", "R", "Y"]:
            parser.error("Wrong list of nucleotides given")

    colors_list = args.color.split(",")
    if len(colors_list) < len(input_files):
        parser.error("Not enough color given")

    if args.legend not in ["True", "False", True, False]:
        print("WARNING : wrong argument given for parameter legend, ")
    if args.legend == "True":
        args.legend = True
    if args.legend == "False":
        args.legend = False

    main(input_files, name_files, nt_list, args.name_fig, args.exon_type, colors_list, args.legend)


if __name__ == "__main__":
    launcher()
