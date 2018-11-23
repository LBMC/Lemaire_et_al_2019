#!/usr/bin/python3

# -*- coding: utf-8 -*-


import numpy
import copy
import sys
import os



def load_fasterDB_coordinates(fasterDB_bed, exon_type):
    """
    Load the coordinates of every fasterDB exon and their surrounding sequences
    :param fasterDB_bed: (string) path to fasterDb bed with 200 nt of intronic sequence before and after the exons
    :param exon_type: (string) the control exon we want to use
    :return: (list of 3 int and  one string) chr, start stop and strand
    """
    list_coordinates = []
    with open(fasterDB_bed, "r") as bed:
        for line in bed:
            if "#" not in line:
                if exon_type == "ALL":
                    line = line.split("\t")
                    list_coordinates.append([line[0], int(line[1]), int(line[2]), line[5]])
                else:
                    if exon_type in line:
                        line = line.split("\t")
                        list_coordinates.append([line[0], int(line[1]), int(line[2]), line[5]])
    return list_coordinates


# def choose_interval(coordinate, size, dic_chr_length):
#     """
#     From a value retrun an interval
#     :param coordinate:  (list of string + 2 int + string) list of coordinates
#     :param size: (int) the size of the interval
#     :param dic_chr_length: (dic of int) length of each chromosomal coordinates
#     :return (list of 2 int) an interval of size ``size``
#     """
#     size = size - 1
#     val = numpy.random.randint(coordinate[1], coordinate[2] + 1)
#     choice = numpy.random.randint(0, 2)
#     if choice == 0:
#         interval =  [val, val + size]
#     else:
#         interval = [val - size, val]
#     count_fail = 0
#     while (interval[0] < 1 or interval[1] > int(dic_chr_length[coordinate[0]])) and count_fail < 10:
#         count_fail += 1
#         choice = numpy.random.randint(0, 2)
#         if choice == 0:
#             interval = [val, val + size]
#         else:
#             interval = [val - size, val]
#     if count_fail == 10:
#         print("WARNING to many failure for interval %s" % str(coordinate))
#     return interval


def choose_interval(coordinate, size, dic_chr_length):
    """
    From a value retrun an interval
    :param coordinate:  (list of string + 2 int + string) list of coordinates
    :param size: (int) the size of the interval
    :param dic_chr_length: (dic of int) length of each chromosomal coordinates
    :return (list of 2 int) an interval of size ``size``
    """
    size = size - 1
    val = numpy.random.randint(coordinate[1], coordinate[2] + 1 - size)
    interval =  [val, val + size]
    return interval


def create_control_coordinates(list_coordinates, nb_control, list_size, dic_chr_length):
    """
    Create control coordinates
    :param list_coordinates: (list of list of string + 2 int + string) list of coordinates
    :param nb_control: (int) the number of control sequence of interest
    :param dic_chr_length :(dic of int) length of each chromosomal coordinates
    :return: (list of string + 2 int + string) list of control coordinates
    """
    print("Create control coordinates ...")
    new_coordinates_list = []
    selected_size = numpy.random.choice(list_size, nb_control, replace=True)
    count = 0
    for i in range(nb_control):
        selected = copy.deepcopy(list_coordinates[numpy.random.randint(0, len(list_coordinates))])
        while (selected[2] - selected[1] + 1) <selected_size[i]:
            selected = copy.deepcopy(list_coordinates[numpy.random.randint(0, len(list_coordinates))])
        prind("Selected coordinates : %s , size %s" % (selected, selected_size[i]))
        new_coord = choose_interval(selected, selected_size[i], dic_chr_length)
        selected[1] = new_coord[0]
        selected[2] = new_coord[1]
        new_coordinates_list.append(selected)
        count += 1
        sys.stdout.write("Coordinates created %s / %s        \r" % (count, nb_control))
    return new_coordinates_list


def create_sequence_dic(hg19_chrom_size_file):
    """
    Create a dictionary from an hg19 fasta files
    :param hg19_chrom_size_file: (string) a file containing the chromosome number and its size
    :return: (dictionary of int) link each chromosome to their size
    """
    dic_records = {}
    with open(hg19_chrom_size_file, "r") as f:
        for line in f:
            line = line.replace("\n", "")
            line = line.split("\t")
            if "_" not in line[0]:
                dic_records[line[0]] = line[1]
    return dic_records


def get_len_list(list_coordinates):
    """
    get the length in nucleotide size of every coordinates in ``list_coordinates``
    :param list_coordinates: (list of string 2 int + str) list of coordinates
    :return: (list of int) list of length of the coordinates in list coordinates
    """
    list_len = []
    count = 0
    for coord in list_coordinates:
        val = int(coord[2]) - int(coord[1]) + 1
        if val > 0:
            list_len.append(val)
        else:
            count += 1
    print("%s sequence with len 0" % count)
    return list_len


def get_sequences_frequencies(fasta_dic, list_size, nb_control):
    print("Calculating control frequencies...")
    list_key = [k for k in list(fasta_dic.keys()) if "GL" not in k]
    ntl = list("ACGT")
    dic_frequencies = {nt: [] for nt in ntl}
    for i in range(nb_control):
        reverse = numpy.random.randint(0, 2)
        chr = numpy.random.choice(list_key)
        size_chosen = numpy.random.choice(list_size)
        max_coord = len(fasta_dic[chr])
        coord1 = numpy.random.randint(0, max_coord)
        dir = numpy.random.randint(0, 2)
        if dir == 0:
            coord2 = coord1 + size_chosen
        else:
            coord2 = coord1 - size_chosen
        while coord2 < 0 or coord2 > max_coord - 1:
            dir = numpy.random.randint(0, 2)
            if dir == 0:
                coord2 = coord1 + size_chosen
            else:
                coord2 = coord1 - size_chosen
        coords = sorted([coord1, coord2])
        if reverse == 1:
            # getting the reverse complement of the sequence
            r = "-"
            sequence = str(fasta_dic[chr][coords[0]:coords[1]].reverse_complement())
        else:
            r = "+"
            sequence = str(fasta_dic[chr][coords[0]:coords[1]])
        while "N" in sequence:
            coord1 = numpy.random.randint(0, max_coord)
            dir = numpy.random.randint(0, 2)
            if dir == 0:
                coord2 = coord1 + size_chosen
            else:
                coord2 = coord1 - size_chosen
            while coord2 < 0 or coord2 > max_coord - 1:
                dir = numpy.random.randint(0, 2)
                if dir == 0:
                    coord2 = coord1 + size_chosen
                else:
                    coord2 = coord1 - size_chosen
            coords = sorted([coord1, coord2])
            if reverse == 1:
                # getting the reverse complement of the sequence
                r = "-"
                sequence = str(fasta_dic[chr][coords[0]:coords[1]].reverse_complement())
            else:
                r = "+"
                sequence = str(fasta_dic[chr][coords[0]:coords[1]])
        prind("sequence of %s, size : %s : %s" % ([chr] + coords + [r], size_chosen, sequence))
        if sequence is not None:
            for nt in ntl:
                dic_frequencies[nt].append((sequence.count(nt) / len(sequence)) * 100 )
        sys.stdout.write("frequencies calculated %s / %s        \r" % (i + 1, nb_control))
    dic_frequencies["S"] = list(numpy.array(dic_frequencies["C"]) + numpy.array(dic_frequencies["G"]))
    dic_frequencies["W"] = list(numpy.array(dic_frequencies["A"]) + numpy.array(dic_frequencies["T"]))
    return dic_frequencies


def get_sequences_frequencies_test(list_coordinates, fasta_dic):
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
        sys.stdout.write("frequencies calculated %s / %s        \r" % (count, total_len))
    dic_frequencies["S"] = list(numpy.array(dic_frequencies["C"]) + numpy.array(dic_frequencies["G"]))
    dic_frequencies["W"] = list(numpy.array(dic_frequencies["A"]) + numpy.array(dic_frequencies["T"]))
    return dic_frequencies


def write_control_file(frequency_dic, filename):
    """
    Write the dictionary ``frequency_dic`` into the file ``filename``
    :param frequency_dic: (dic of list of float) a dictionary containing the frequencies of many control sequences
    :param filename: (string) file where the control dic will be written
    """
    with open(filename, "w") as f:
        f.write("ctrl_dic = %s\n" % str(frequency_dic))


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


def control_dictionaries_generator(control_sequence_dic, nb_iteration, size):
    """
    Generate a mean dictionary containing for each nucleotide their ``nb_iteration`` mean frequencies \
     calculated using a sample of ``size`` control sequence
    :param control_sequence_dic: (dic of list of float) for each iupac nucleotides, give the frequencies of the \
    control sequences selected
    :param nb_iteration: (int) the number of time we gonna select ``size`` frequencies of each nt
    :param size: (int) the size of each control set created
    :return: (dic of list of float) each list of the dictionary corresponds to the mean frequencies of a nt \
     for each control set selected
    """
    print("Creating control sample...")
    max_val = len(control_sequence_dic["A"])
    dic_mean_control_set = {nt: [] for nt in control_sequence_dic.keys()}
    for i in range(nb_iteration):
        control_sample =  {nt: [] for nt in control_sequence_dic.keys()}
        for j in range(size):
            val = numpy.random.randint(0, max_val)
            for nt in control_sample.keys():
                control_sample[nt].append(control_sequence_dic[nt][val])

        prind("Sample %s : %s" % (i, str(control_sample)))
        for nt in dic_mean_control_set.keys():
            dic_mean_control_set[nt].append(numpy.mean(control_sample[nt]))
        sys.stdout.write("Median sample done : %s / %s             \r" % (i, nb_iteration))
    prind("result mean : %s" % str(dic_mean_control_set))
    return dic_mean_control_set


def make_control_sets(nb_iteration, hg19_chromosome_file, fasterdb_bed, exon_type, test_coordinates, fasta_dic):
    output = os.path.realpath(os.path.dirname(__file__)) + "/control/"
    if not os.path.isdir(output):
        os.mkdir(output)
    output_file = "%s%s_control_frequencies.py" % (output, exon_type)
    if not os.path.isfile(output_file):
        list_size = get_len_list(test_coordinates)
        nb_control = 1000000
        # nb_control = 5
        dic_frequency = get_sequences_frequencies(fasta_dic, list_size, nb_control)
        prind("frequencies that will be written : \n ----> %s " % str(dic_frequency))
        write_control_file(dic_frequency, output_file)
    else:
        sys.path.insert(0, output)
        mod = __import__(os.path.basename(output_file).replace(".py", ""))
        dic_frequency = mod.ctrl_dic
    dic_frequencies = control_dictionaries_generator(dic_frequency, nb_iteration, len(test_coordinates))
    # dic_frequencies = control_dictionaries_generator(dic_frequency, 100, 4)
    # dic_frequencies = dic_frequency
    return dic_frequencies


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