#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    Launches the ``launcher script`` for every file in the Rename file
"""

import subprocess
import os
import argparse

def find_files(input_bed):
    """
    Return the files name contained in ``input_bed``
    :param input_bed:  (string) folder containing bed files
    :return: (list of string) files contained in ``input_bed`` folder
    """
    cmd = "ls %s" % input_bed
    list_files = subprocess.check_output(cmd.split(" "), stderr=subprocess.STDOUT).decode("ascii")
    list_files = list_files.split("\n")
    list_files = list_files[0:-1]
    new_list = []
    for nfile in list_files:
        new_list.append(input_bed + nfile)
    return new_list[0:-1]


def launcher_program(input_bed, output_folder):
    """
    Launch the launcher program
    :param input_bed: (string) bed file
    :param output_folder: (string) output where the result will be created
    """
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    cmd =  "python3 src/launcher.py --clip_bed %s --output %s --nb_iteration 1000 --enrichment Y" % (input_bed,
                                                                                                     output_folder)
    print(cmd)
    subprocess.check_call(cmd.split(" "), stderr=subprocess.STDOUT)


def main(bed_folder, output_folder):
    """
    Launches the launcher.py file for every bed in a folder
    :param bed_folder: (string) path to bed folder
    :param output_folder: (string) file where the result will be created
    """
    list_files = find_files(bed_folder)
    for my_file in list_files:
        base_name = os.path.basename(my_file).split(".")[0]
        output = output_folder + base_name + '/'
        if not os.path.isfile("%sfigures/A_density_%s.html" % (output, base_name)):
            launcher_program(my_file, output)
        else:
            print("File %s already analyzed" % my_file)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    This program aims to launch the ``launcher.py`` for every bed files located in a folder given in input
    """)
    # Arguments for the parser

    req_arg = parser.add_argument_group("Required arguments")
    req_arg.add_argument("--input", dest="input", help="a folder containing bed files", required=True)

    req_arg.add_argument("--output", dest="output", help="a folder containing where the result will be created files", required=True)
    args = parser.parse_args()

    if not os.path.isdir(args.output):
        parser.error("The output directory doesn't exist !")
    if args.output[-1] != "/":
        args.output += "/"

    if not os.path.isdir(args.input):
        parser.error("The input directory doesn't exist !")
    if args.input[-1] != "/":
        args.input += "/"

    main(args.input, args.output)


if __name__ == "__main__":
    launcher()
