#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import subprocess
import gzip
import os


def get_files(folder):
    """
    List the files in ``folder``
    :param folder: (string) a folder
    :return: (list of string) the files located in ``folder``
    """
    files = os.listdir(folder)
    return [folder + namefile for namefile in files]


def bed_0based_2_1based(bed_file, output):
    """
    Convert a 0 based bed into  a 1 based bed
    :param bed_file: (string) a bed obtain from a clip experiment
    :param output: (string) path where the 1 based bed will be written
    :return: (string) the name of the
    """
    basename = bed_file.split("/")[-1].replace(".gz", "").replace(".bed", "")
    new_name = "%s%s.1b.bed.gz" % (output, basename)
    if "gz" in bed_file:
        bedin = gzip.open(bed_file, "rt")
    else:
        bedin = open(bed_file, "r")
    with gzip.open(new_name, "wt") as bedout:
        for line in bedin:
            if line[0] == "#":
                bedout.write(line)
            else:
                line = line.split("\t")
                line[0] = line[0].replace("chr", "")
                if line[0] == "M":
                    line[0] += "T"
                line[1] = str(int(line[1]) + 1)
                bedout.write("\t".join(line))
    bedin.close()
    return new_name


def intersect_bed(faster_db_bed, clib_1b_bed, output, name):
    """
    Keep only the feature of ``clib_1b_bed``` if they overlap a feature in ``faster_db_bed``
    :param faster_db_bed: (string) a bed containing every fasterDB exons with 200bp of interval
    :param clib_1b_bed:(string) a bed obtain from a clip experiment 1 based
    :param output: (string) path where the intersect bed will be written
    :param name: (string) the name of the file
    :return: (string) the name of the intersect bed
    """
    basename = clib_1b_bed.split("/")[-1].replace(".1b.bed.gz", "")
    new_name = "%s%s_Intersect_%s.1b.bed.gz" % (output, basename, name)
    cmd = "intersectBed -a %s -b %s -u -f 1 | gzip -c > %s " % (clib_1b_bed, faster_db_bed, new_name)
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
    return new_name


def meme_launcher(meme_path, fasta_file, output):
    """
    Create file that will be use to create the enriched motifs

    :param meme_path: (string) path to the folder containing the meme programs
    :param fasta_file: (string) path to the fasta file containing the clipped sequence
    :param output: (string) path where the result will be created
    """
    cmd1 = "%sfasta-center -dna -len 100 < %s 1> %sseqs-centered.fa" % (meme_path, fasta_file, output)
    cmd2 = "%sfasta-get-markov -nostatus -nosummary -dna -m 1 %s %sbackground" % (meme_path, fasta_file, output)
    cmd3 = "%smeme %sseqs-centered.fa -oc %smeme_weblogo -mod zoops -nmotifs 1 -minw 6 -maxw 30 -bfile %sbackground " \
           "-dna -revcomp -nostatus" % (meme_path, output, output, output)
    # cmd4 = "%sdreme -verbosity 1 -oc %sdreme_weblogo -png -dna -p %sseqs-centered.fa -m 3 -k 9 -norc" % \
    #        (meme_path, output, output)
    print(cmd1)
    subprocess.call(cmd1, shell=True, stderr=subprocess.STDOUT)
    print(cmd2)
    subprocess.call(cmd2, shell=True, stderr=subprocess.STDOUT)
    print(cmd3)
    subprocess.call(cmd3, shell=True, stderr=subprocess.STDOUT)
    # print(cmd4)
    # subprocess.call(cmd4, shell=True, stderr=subprocess.STDOUT)


def get_zagros_weblogo(zagros_file, output):
    """
    Creation of a weblogo thanks to zagros results
    :param zagros_file: (string) file containing zagros result
    :param output: (string) ath where the weblogo will be created
    """
    list_motif = []
    motif = False
    list_seq = None
    with open(zagros_file, "r") as zf:
        for line in zf:
            line = line.replace("\n", "").split("\t")
            if not motif:
                if line[0] == "BS":
                    motif = True
                    list_seq = [line[1].split(";")[0]]
            else:
                if line[0] == "BS":
                    list_seq.append(line[1].split(";")[0])
                else:
                    list_motif.append(list_seq)
                    del list_seq
                    motif = False
    mod = __import__("figure_maker")
    for i in range(len(list_motif)):
        mod.web_logo_creator(list_motif[i], "zagros_motif_%s" % i, output)


def zagros_launcher(zagros_path, fasta_file, size, output):
    """
    Launch zagros program to predict rna binding sites.

    :param zagros_path: (string) path where the zagros bins are located
    :param fasta_file: (string) path where the peak are located
    :param size: (int) the size of the motif that we want zagros to find
    :param output: (string) path where the zagros result will be created
    """
    output += "zagros/"
    if not os.path.isdir(output):
        os.mkdir(output)
    res_file = "%szagros_result.txt" % output
    cmd1 = "%sthermo -o %sinput.str %s" % (zagros_path, output, fasta_file)
    cmd2 = "%szagros -t %sinput.str -o %s -w %s -n 3 %s" % (zagros_path, output, res_file, size, fasta_file)
    print(cmd1)
    subprocess.call(cmd1, shell=True, stderr=subprocess.STDOUT)
    print(cmd2)
    subprocess.call(cmd2, shell=True, stderr=subprocess.STDOUT)
    get_zagros_weblogo(res_file, output)
