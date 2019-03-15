#!/usr/bin/env python3

"""
Description :
    The goal of this script is to compute the number of stretches for each sequences located before an exon.
"""

import sys
import config


def stretch_finder_nt(sequence, nt, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param nt: (string) the name of the nt to return
    :param stretch_len: (int) the length of the stretch of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: the number of stretch of the nucleotide "nt" here
    """
    nb_stretch = 0
    if nt in ["A", "T", "G", "C"]:
        for i in range(len(sequence) - stretch_len + 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter == nt:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    else:
        iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
        for i in range(len(sequence) - stretch_len + 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter in iupac[nt]:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    return nb_stretch


def stretch_counter(exon_list, stretch_data, sequence_boundaries):
    """
    Launch stretch_finder_nt  for every exon in the exon list.
    :param exon_list: (list of ExonClass object) list of exons
    :param stretch_data: (list of 2 int) list containing the total length and the size content of the strectch
    :param sequence_boundaries: (string) the indice of the upstream sequence wanted
    :return: (dictionary of list of int) each nucleotide (key of the dictionary) is linked to
    """
    stretch_dic = {nt: [] for nt in config.nt_list}
    my_len = len(exon_list)
    for nt in stretch_dic.keys():
        count = 0
        for exon in exon_list:
            count += 1
            if exon.upstream_intron.sequence_proxi is not None:
                if count != my_len:
                    sys.stdout.write("\t\t--> recovering %s stretches (%s/%s) : %s / %s                  \r" %
                                     (nt, stretch_data[1], stretch_data[0], count, my_len))
                else:
                    sys.stdout.write("\t\t--> recovering %s stretches (%s/%s) : %s / %s                  \n" %
                                     (nt, stretch_data[1], stretch_data[0], count, my_len))
                sys.stdout.flush()
                sequence = exon.upstream_intron.sequence_proxi[sequence_boundaries[0]:sequence_boundaries[1]]
                stretch_dic[nt].append(stretch_finder_nt(sequence, nt, stretch_data[0], stretch_data[1]))
    return stretch_dic
