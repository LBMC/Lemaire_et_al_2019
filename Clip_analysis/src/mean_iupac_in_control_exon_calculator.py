#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    This script aims to calculates for a particular exon type the mean frequencies of very nucleotides in this set.
    Tis also aims to calculates the mean iupac content of a particular bed file
"""

import sqlite3
import numpy
import get_bed_sequences
import os
import sys

def get_iupac_exon(cnx, exon_type):
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = "SELECT iupac_exon FROM sed WHERE exon_type LIKE '%{}%'".format(exon_type)
    else:
        query = "SELECT iupac_exon FROM sed"
    cursor.execute(query)
    res = cursor.fetchall()
    dic_res = {nt: [] for nt in list("ACGTSW")}
    for iupac in res:
        if iupac[0]:
            dic_res["A"].append(float(iupac[0].split(";")[0]))
            dic_res["C"].append(float(iupac[0].split(";")[1]))
            dic_res["G"].append(float(iupac[0].split(";")[2]))
            dic_res["T"].append(float(iupac[0].split(";")[3]))
            dic_res["S"].append(float(iupac[0].split(";")[4]))
            dic_res["W"].append(float(iupac[0].split(";")[5]))
    return dic_res


def get_sequences_frequencies(list_coordinates, fasta_dic):
    """
    Get the nucleotide frequency in sequences defined by each set of coordinates in ``list_coordinates``
    :param list_coordinates: (list of 3 ints) the chromosome and the chromosomal coordinates of each features
    :param fasta_dic: (dictionary of Seq object) dictionary where the chromomsome number are the key and their sequence\
     are defined in the values of the dictionary
    :return: (dic of list of float) for each iupac nucleotides, give the frequencies of the sequences defined by \
    ``list_coordinates``
    """
    print("Calculating control frequencies...")
    iupac = list("ACGTSW")
    ntl = list("ACGT")
    dic_frequencies = {nt: [] for nt in ntl}
    count = 0
    total_len = len(list_coordinates)
    for coord in list_coordinates:
        if coord[-1] == "-":
            # getting the reverse complement of the sequence
            sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]].reverse_complement())
        else:
            sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]])
        prind("sequence of %s : %s" % (coord, sequence))
        if sequence is not None:
            for nt in ntl:
                dic_frequencies[nt].append((sequence.count(nt) / len(sequence)) * 100 )
            # dic_frequencies["S"] = dic_frequencies["C"] + dic_frequencies["G"]
            # dic_frequencies["W"] = dic_frequencies["A"] + dic_frequencies["T"]
        count += 1
    dic_frequencies["S"] = list(numpy.array(dic_frequencies["C"]) + numpy.array(dic_frequencies["G"]))
    dic_frequencies["W"] = list(numpy.array(dic_frequencies["A"]) + numpy.array(dic_frequencies["T"]))
    return dic_frequencies


def create_mean_frequency_dic(frequencies_dic):
    """
    From a dictionary of list of float return a mean dictionary for each keys
    :param frequencies_dic: (dic of list of float) each key of the dictionary corresponds to an nucleotide and its \
    linked to their frequencies in control or test sequence
    :return: (dic of float) mean frequencies of each nucleotide
    """
    meean_dic = {}
    for nt in frequencies_dic.keys():
        meean_dic[nt] = numpy.mean(frequencies_dic[nt])
    prind("Mean dic : %s" % str(meean_dic))
    return meean_dic


def set_debug(debug=0):
    """
    Set the debug mode
    :param debug: (int) 0 debug mode disabled, 1 enabled
    """
    global d
    d = debug


def prind(message):
    """
    print a message in debug mode
    :param message: (string) the message to print
    """
    if d == 1:
        print(message)


def main():
    exon_type = "ALL"
    set_debug(0)
    hg19_fasta = os.path.realpath(os.path.dirname(__file__)).replace("src","data/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src","data/sed.db")
    cnx = sqlite3.connect(seddb)
    bed_file = sys.argv[1]
    exon_dic = get_iupac_exon(cnx, exon_type)
    ctrl_dic = create_mean_frequency_dic(exon_dic)
    print("Iupac frequencies for %s exons : " % exon_type)
    print(str(ctrl_dic))
    list_coordinates = get_bed_sequences.get_bed_coordinates(bed_file)
    print("loading hg19 data !")
    dic_records = get_bed_sequences.create_sequence_dic(hg19_fasta)
    freq_bed = get_sequences_frequencies(list_coordinates, dic_records)
    final_bed = create_mean_frequency_dic(freq_bed)
    print("Mean frequencies in  bed :")
    print(str(final_bed))


if __name__ == "__main__":
    main()