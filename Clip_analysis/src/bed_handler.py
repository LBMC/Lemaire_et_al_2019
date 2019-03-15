#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import subprocess
import gzip
import os
import union_dataset_function
import fasterdb_bed_add_exon_type


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


def bed_creator(cnx, cnx_sed, dest_folder, sf_name, regulation, chrom_size_file):
    """
    Create a bed file containing all the exons regulated by ``sf_name`` with ``regulation``.

    :param cnx: (sqlite3 connector object) connection to fasterDB lite database
    :param cnx: (sqlite3 connector object) connection to sed database
    :param dest_folder: (string) path where the bed will be created
    :param sf_name: (string) the name of a splicing factor
    :param regulation: (string) up or down
    :param: (string) a file containing chromosome size
    :return: (string) the name of the bed file created
    """
    sf_name = sf_name.upper()
    sf_name = sf_name.replace("SFRS", "SRSF")
    dic_fact = {"TRA2A": "TRA2A_B"}
    strand_dic = {-1:"-", 1:"+"}
    if sf_name in dic_fact:
        exon_list = union_dataset_function.get_every_events_4_a_sl(cnx_sed, dic_fact[sf_name], regulation)
    else:
        exon_list = union_dataset_function.get_every_events_4_a_sl(cnx_sed, sf_name, regulation)
    if len(exon_list) == 0:
        print("%s exon %s : %s %s" % ("\033[0;31m", sf_name, len(exon_list), "\033[0m"))
    else:
        print("%s exon %s : %s %s" % ("\033[0;32m", sf_name, len(exon_list), "\033[0m"))
    cursor = cnx.cursor()
    exon_info = []
    for exon in exon_list:
        query = """SELECT t1.chromosome, t1.start_on_chromosome, t1.end_on_chromosome, t2.official_symbol, t1.pos_on_gene,
                   t2.strand, t1.end_on_chromosome - t1.start_on_chromosome + 1
                   FROM exons t1, genes t2
                   WHERE t1.id_gene = t2.id
                   AND t1.id_gene = %s
                   AND t1.pos_on_gene = %s
                   ORDER BY t1.chromosome ASC, t1.start_on_chromosome ASC
                """ % (exon[0], exon[1])
        cursor.execute(query)
        res = cursor.fetchall()
        if len(res) > 1:
            print("Error, only one exon should be found for %s_%s exon" %  (exon[0], exon[1]))
            exit(1)
        my_exon = list(res[0][0:3]) + ["%s_%s" % (res[0][3], res[0][4])] + ["."] + [strand_dic[res[0][5]], res[0][6]]
        exon_info.append("\t".join(list(map(str, my_exon))))
    bed_content = "\n".join(exon_info)
    filename = "%sSF_%s_exons-union.bed" % (dest_folder, regulation)
    final_name = filename.replace(".bed", "_add200nt.bed")
    with open(filename, "w") as bedfile:
        bedfile.write(bed_content + "\n")
    fasterdb_bed_add_exon_type.add_intron_sequence(filename, final_name, chrom_size_file)
    return final_name


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
