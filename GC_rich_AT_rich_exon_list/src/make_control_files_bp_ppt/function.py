#!/sur/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    From a list of exon containing the upstream proxi sequence (100 nt before the exon, find the best \
    branch point and ppt score
"""


import subprocess
import pandas as pd
from io import StringIO
import os
import sys

# global variable
svm_launcher = "/media/nicolas/DD_1/Splicing_Lore_project/SVM_BP_finder/svm_bpfinder.py"
output = os.path.realpath(os.path.dirname(__file__)).replace("src/make_control_files_bp_ppt", "result/bp_ppt_score/")


def fasta_writer(exon, output, size):
    """
    Write a file named input.fa in the directory ``output``
    :param exon: (ExonClass exon) an exons
    :param output: (string) path where the output file will be created
    :param size: (string) the size of the intron wanted
    :return: (string) the name of the input file
    """
    if size == 100:
        sequence = exon.upstream_intron.sequence_proxi
    else:
        sequence = exon.upstream_intron.sequence_proxi[-size:]
    input_file = "%sinput.fa" % output
    with open(input_file, "w") as out_file:
        out_file.write(">%s_%s_proxi_intron_seq\n" % (exon.gene.name, exon.position))
        out_file.write(sequence)
    return(input_file)


def svm_bp_finder_launcher(svm_bp_launcher, input_file):
    """
    Launch svm_bp_finder
    :param svm_bp_launcher: (string)svm_bp finder launcher
    :param input_file: (string) input fasta file
    :return: (2 floats) branch point score and ppt score
    """
    result = subprocess.check_output([svm_bp_launcher, "--input", input_file,
                                      "--species", "Hsap", "-l", "100"]).decode("ascii")
    result = StringIO(result)
    df = pd.read_csv(result, sep="\t")
    subprocess.check_call(['rm', input_file])
    if not df.empty:
        nb_putative_bp = df.shape[0]
        nb_good_bp = sum(df["svm_scr"] > 0)
        df = df.sort_values(by=["svm_scr"], ascending=False).head(1)
        bp_score = float(df["bp_scr"])
        ppt_score = float(df["ppt_scr"])
        return bp_score, ppt_score, nb_putative_bp, nb_good_bp
    else:
        return None, None, 0, 0


def bp_ppt_calculator(exon_list, size=100):
    """
    Launch svm_bp finder for every exon in the exon list
    :param exon_list: (list of ExonClass object) list of exons
    :param size: (string) the size of the upstream sequence wanted
    :return: (2 lists of floats) lists of float of bp
    """
    bp_score_list = []
    ppt_score_list = []
    nb_bp_list = []
    nb_good_bp_list = []
    if not os.path.isdir(output):
        os.mkdir(output)
    count = 0
    exon_len = len(exon_list)
    for exon in exon_list:
        count += 1
        sys.stdout.write("%s / %s      \r" % (count, exon_len))
        if exon.upstream_intron.sequence_proxi is not None:
            input_file = fasta_writer(exon, output, size)
            bp_score, ppt_score, nb_bp, nb_good_bp = svm_bp_finder_launcher(svm_launcher, input_file)
            if bp_score is not None:
                bp_score_list.append(bp_score)
            if ppt_score is not None:
                ppt_score_list.append(ppt_score)
            if nb_bp is not None:
                nb_bp_list.append(nb_bp)
            if nb_good_bp is not None:
                nb_good_bp_list.append(nb_good_bp)
    return bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list


