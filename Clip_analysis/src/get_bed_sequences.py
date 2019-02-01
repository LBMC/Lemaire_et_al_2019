#!/urs/bin/env python3

# -*- coding utf-8 -*-

"""
Description:

    From a given bed file return the wanted sequence of each element in the bed file
"""

import gzip
from Bio import SeqIO
import math


# Functions

def read_gzip_bed_file(bed_file, list_chromosome):
    """
    Get the content of a compressed bed file
    :param bed_file: (string) the name of a gunzipped bed file
    :param list_chromosome: (list of string) list of chromosome authorized to retrieve coordinates
    :return: (list of list of 3 ints and 1 string) the chromosome and the chromosomal coordinates and \
     the strand of each features
    """
    list_coordinates = []
    with gzip.open(bed_file, "rt") as bed:
        for line in bed:
            if line[0] != "#":
                line = line.split("\t")
                if line[0] in list_chromosome:
                    list_coordinates.append([line[0], int(line[1]), int(line[2]), line[5]])
    return list_coordinates


def read_bed_file(bed_file, list_chromosome):
    """
    Get the content of a compressed bed file
    :param bed_file: (string) the name of a bed file
    :param list_chromosome: (list of string) list of chromosome authorized to retrieve coordinates
    :return: (list of list of 3 ints and 1 string) the chromosome and the chromosomal coordinates and \
     the strand of each features
    """
    list_coordinates = []
    with open(bed_file, "r") as bed:
        for line in bed:
            if "#" not in line:
                line = line.split("\t")
                if line[0] in list_chromosome:
                    list_coordinates.append([line[0], int(line[1]), int(line[2]), line[5]])
    return list_coordinates


def get_bed_coordinates(bed_file, list_chromosome):
    """
    Get the coordinates of a bed files
    :param bed_file: (string) a bed files (that can be compressed or not)
    :param list_chromosome: (list of string) list of chromosome authorized to retrieve coordinates
    :return: (list of list of 3 ints and 1 strings) the chromosome and the chromosomal coordinates and \
     the strand of each features
    """
    if ".gz" in bed_file:
        list_coordinates = read_gzip_bed_file(bed_file, list_chromosome)
    else:
        list_coordinates = read_bed_file(bed_file, list_chromosome)
    return list_coordinates


def clip_sequence_size(list_coordinates):
    """
    Return the range of the len of clip peaks
    :param list_coordinates:  the chromosome and the chromosomal coordinates and \
     the strand of each features
    :return: (2 int) the max and the min size of the clips
    """
    min_val = list_coordinates[0][1] - list_coordinates[0][2] + 1
    max_val = min_val
    for coord in list_coordinates:
        len_seq = coord[2] - coord[1] + 1
        if len_seq < min_val:
            min_val = len_seq
        if len_seq > max_val:
            max_val = len_seq
    return min_val, max_val


def get_sequences(list_coordinates, fasta_dic):
    """
    Get the sequence defined by each set of coordinates in ``list_coordinates``
    :param list_coordinates: (list of 3 ints) the chromosome and the chromosomal coordinates of each features
    :param fasta_dic: (dictionary of Seq object) dictionary where the chromomsome number are the key and their sequence\
     are defined in the values of the dictionary
    :return: (list of strings) list of sequence of interest
    """
    list_seq = []
    for coord in list_coordinates:
        if coord[-1] == "-":
            # getting the reverse complement of the sequence
            sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]].reverse_complement())
        else:
            sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]])
        if len(sequence) > 5:
            list_seq.append(sequence)
    return list_seq


def totaly_defined(sequence):
    """
    Says if all the nucleotide are well defined within the sequence

    :param sequence: (string) nucleotide sequence
    :return: (boolean) True if all nucleotide are well defined, False else
    """
    seq_defined = sequence.count("A") + sequence.count("T") + sequence.count("G") + sequence.count("C")
    if seq_defined == len(sequence):
        return True
    return False


def get_list_sequence_hexanucleotide_frequencies(list_seq):
    """
    Get the average frequencies of every hexanucleotide in a list of nucleotide sequence ``list_seq``
    :param list_seq: (list of string) list of nucleotide sequence
    :return: (dictionary of float) dictionary that containing the mean value for each hexanucleotide in ``list_seq``
    """
    mean_hexa_count = {}
    nb_sequence_selected = 0
    for sequence in list_seq:
        if totaly_defined(sequence):
            hex_dic = get_hexanucleotide_frequencies(sequence)
            for hexa in hex_dic.keys():
                if hexa not in mean_hexa_count.keys():
                    mean_hexa_count[hexa] = [hex_dic[hexa]]
                else:
                    mean_hexa_count[hexa].append(hex_dic[hexa])
            nb_sequence_selected += 1
    final_dic = {}
    for hexa in mean_hexa_count.keys():
        final_dic[hexa] = sum(mean_hexa_count[hexa]) / nb_sequence_selected
    return final_dic


def get_hexanucleotide_frequencies(sequence):
    """
    Get the frequencies of each hexanucleotide within the nucleotide sequence ``sequence``
    :param sequence: (string) a nucleotide sequence
    :return: (dictionary of float) link each hexanucleotide to it's frequency within ``sequence``
    """
    dic_freq = {}
    for i in range(len(sequence) - 5):
        hexa = sequence[i:i+6]
        if hexa not in dic_freq:
            dic_freq[hexa] = 1 / (len(sequence) - 5)
        else:
            dic_freq[hexa] += 1 / (len(sequence) - 5)
    return dic_freq


def write_fasta_file(list_coordinates, fasta_dic, output):
    """
    Cerate a fasta file that will be used by meme to find an enriched motif in those sequence
    :param list_coordinates: (list of 3 ints + ) the chromosome and the chromosomal coordinates of each features
    :param fasta_dic:(dictionary of Seq object) dictionary where the chromomsome number are the key and their sequence\
     are defined in the values of the dictionary
    :param output: (string) path where the fasta file will be created
    """
    name_file = "%sclip_peaks.fa" % output
    with open(name_file, "w") as fasta:
        count = 0
        for coord in list_coordinates:
            count += 1
            if coord[-1] == "-":
                # getting the reverse complement of the sequence
                sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]].reverse_complement())
            else:
                sequence = str(fasta_dic[coord[0]][coord[1] - 1:coord[2]])
            if len(sequence) > 5:
                fasta.write(">seq%s|1\n" % count)
                fasta.write(sequence + "\n")
    return name_file


def create_sequence_dic(hg19_fasta):
    """
    Create a dictionary from an hg19 fasta files
    :param hg19_fasta: (string) a fasta files containing the sequence of every chromosome in hg19
    :return: (dictionary of Seq object) dictionary where the chromomsome number are the key and their sequence\
     are defined in the values of the dictionary
    """
    dic_records = {}
    for record in SeqIO.parse(hg19_fasta, "fasta"):
        dic_records[record.id] = record.seq
    return dic_records
