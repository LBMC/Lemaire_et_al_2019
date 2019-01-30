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
sequence_u2 = "GTGTAGTA"
bound = {"A": 2, "T": 2, "G": 3, "C": 3}


def fasta_writer(exon, output_folder, size):
    """
    Write a file named input.fa in the directory ``output``
    :param exon: (ExonClass exon) an exons
    :param output_folder: (string) path where the output file will be created
    :param size: (string) the size of the intron wanted
    :return: (string) the name of the input file
    """
    if size == 100:
        sequence = exon.upstream_intron.sequence_proxi
    else:
        sequence = exon.upstream_intron.sequence_proxi[-size:]
    input_file = "%sinput_bp_finder.fa" % output_folder
    with open(input_file, "w") as out_file:
        out_file.write(">%s_%s_proxi_intron_seq\n" % (exon.gene.name, exon.position))
        out_file.write(sequence)
    return input_file


def adapt_sequence(sequence, align):
    """
    Remove every position with a dot in align from sequence.

    :param sequence: (string) a nucleotide sequence
    :param align: (string) list of ( or ) and dot
    :return: (string) the sequence without any position corresponding to a dot in ``align``
    """
    s = list(sequence)
    adjust = 0
    for i in range(len(align)):
        if align[i] == ".":
            del s[i - adjust]
            adjust += 1
    return "".join(s)


def rna_duplex_launcher(sequence):
    """
    Launch RNA duplex program.

    :return: (int) The number of hydrogen bound in the hybridization of sequence and sequence_u2
    """
    cmd = """RNAduplex <<EOF
>seq1
%s
>seq2
%s
EOF""" % (sequence, sequence_u2)
    res = subprocess.check_output(cmd, shell=True).decode("ascii")
    pos_seq1 = None
    pos_seq2 = None
    line = res.split("\n")[-2]
    line = list(filter(None, line.split(" ")))
    try:
        pos_seq1 = list(map(int, line[1].split(",")))
        pos_seq2 = list(map(int, line[3].split(",")))
    except ValueError:
        print("ERROR : wrong position values : %s " % line[1])
        print("ERROR : wrong position values : %s " % line[3])
        exit(1)
    align = line[0].split("&")
    if align[0].count(".") + align[0].count("(") != len(align[0]):
        print("Error : not only ( and . symbol in align : %s" % align[0])
        exit(1)
    if align[1].count(".") + align[1].count(")") != len(align[1]):
        print("Error : not only ) and . symbol in align : %s" % align[1])
        exit(1)
    if pos_seq1[1] == pos_seq1[0]:
        return 0

    new_seq1 = sequence[pos_seq1[0]-1:pos_seq1[1]]
    new_seq1 = adapt_sequence(new_seq1, align[0])
    new_seq2 = sequence_u2[pos_seq2[0]-1:pos_seq2[1]]
    new_seq2 = adapt_sequence(new_seq2, align[1])
    new_seq2 = new_seq2[::-1]
    if len(new_seq1) != len(new_seq2):
        print("Error, the sequence length differ : s1 : %s, s2: %s" % (new_seq1, new_seq2))
    hbound = 0
    for i in range(len(new_seq1)):
        hbound += min(bound[new_seq1[i]], bound[new_seq2[i]])
    return hbound


def motif_calculator(sequence, motif):
    """
    Calculate the frequencies of AG di-nucleotide downstream of the branch point.

    :param sequence: (string) sequence downstream of the branch point
    :param motif: (string) the motif we want to count
    :return: (float) the AG dinucleotide frequencies
    """
    if len(sequence) >= len(motif):
            return sequence.count(motif)
    return None


def svm_bp_finder_launcher(svm_bp_launcher, input_file, exon, size):
    """
    Launch svm_bp_finder
    :param svm_bp_launcher: (string)svm_bp finder launcher
    :param exon: (ExonClass instance) an exon
    :param input_file: (string) input fasta file
    :param size: (string) the size of the sequence
    :return: (2 floats, 2 int, and 1 string) branch point score and ppt score, nb_putative branch point and sequence
    """
    uaa_motif = motif_calculator(exon.upstream_intron.sequence_proxi[-size:], "TAA")
    una_motif = 0
    for motif in ["TAA", "TCA", "TGA", "TTA"]:
        una_motif += motif_calculator(exon.upstream_intron.sequence_proxi[-size:], motif)
    result = subprocess.check_output([svm_bp_launcher, "--input", input_file,
                                      "--species", "Hsap", "-l", "100"]).decode("ascii")
    result = StringIO(result)
    df = pd.read_csv(result, sep="\t")
    subprocess.check_call(['rm', input_file])
    if not df.empty:
        nb_putative_bp = df.shape[0]
        nb_good_bp = sum(df["svm_scr"] > 0)
        df = df.sort_values(by=["svm_scr"], ascending=False).head(1)
        if float(df["svm_scr"]) > 0:
            seq = exon.upstream_intron.sequence_proxi[-size:]
            bp = int(df["ss_dist"])
            ag_count = motif_calculator(seq[-bp + 1:], "AG")
            if len(seq) >= bp + 3 and bp > 8:
                start = int(-bp-3)
                stop = int(-bp+9)
                sequence = seq[start:stop]
            else:
                print("Warning sequence length to short: bp pos %s, sequence : %s, len : %s" % (bp, seq, len(seq)))
                sequence = None
            if bp > 3 and len(seq) >= bp + 5:
                start = int(-bp-5)
                stop = int(-bp + 4)
                sequ2 = seq[start:-bp] + seq[-bp+1:stop]
                hbound = rna_duplex_launcher(sequ2)
            else:
                hbound = None
        else:
            sequence = None
            ag_count = None
            hbound = None
        bp_score = float(df["bp_scr"])
        ppt_score = float(df["ppt_scr"])
        return bp_score, ppt_score, nb_putative_bp, nb_good_bp, sequence, ag_count, hbound, uaa_motif, una_motif
    else:
        return None, None, 0, 0, None, None, None, uaa_motif, una_motif


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
    sequence_list = []
    ag_count_list = []
    hbound_list = []
    uaa_motif_list = []
    una_motif_list = []
    if not os.path.isdir(output):
        os.mkdir(output)
    count = 0
    exon_len = len(exon_list)
    for exon in exon_list:
        count += 1
        sys.stdout.write("%s / %s      \r" % (count, exon_len))
        if exon.upstream_intron.sequence_proxi is not None:
            input_file = fasta_writer(exon, output, size)
            bp_score, ppt_score, nb_bp, nb_good_bp, sequence, ag_count, hbound, \
                uaa_motif, una_motif = svm_bp_finder_launcher(svm_launcher, input_file, exon, size)
            if bp_score is not None:
                bp_score_list.append(bp_score)
            if ppt_score is not None:
                ppt_score_list.append(ppt_score)
            if nb_bp is not None:
                nb_bp_list.append(nb_bp)
            if nb_good_bp is not None:
                nb_good_bp_list.append(nb_good_bp)
            if sequence is not None:
                sequence_list.append(sequence)
            if ag_count is not None:
                ag_count_list.append(ag_count)
            if hbound is not None:
                hbound_list.append(hbound)
            if uaa_motif is not None:
                uaa_motif_list.append(uaa_motif)
            if una_motif is not None:
                una_motif_list.append(una_motif)
    return bp_score_list, ppt_score_list, nb_bp_list, nb_good_bp_list, sequence_list, \
        ag_count_list, hbound_list, uaa_motif_list, una_motif_list
