#!/sur/bin/env python3

# -*- coding: utf-8 -*-

"""
Description:
    From a list of exon containing the upstream proxi sequence (100 nt before the exon, find the best \
    branch point and ppt score
"""


import subprocess
import os


def fasta_writer(exon_list, filename, splicing_site):
    """
    Write a file named input.fa in the directory ``output``
    :param exon_list: ( list of ExonClass) list of exons
    :param filename: (string) list of exons
    :param splicing_site: (string) 3ss or 5ss
    :return: (string) the name of the input file
    """
    count = 0
    with open(filename, "w") as out_file:
        for exon in exon_list:
            if splicing_site == "3ss":
                if exon.seq_3ss is not None:
                    out_file.write(">seq_%s_%s\n" % (exon.gene.name.replace(".", ","), exon.position))
                    out_file.write(exon.seq_3ss + "\n")
                else:
                    print("WARNING : exons %s_%s with %s start_on_gene can't be written : gene_len : %s" %
                          (exon.gene.name, exon.position, exon.start, exon.gene.length))
                    count += 1
            else:
                if exon.seq_5ss is not None:
                    out_file.write(">seq_%s_%s\n" % (exon.gene.name.replace(".", ","), exon.position))
                    out_file.write(exon.seq_5ss + "\n")
                else:
                    print("WARNING : exons %s_%s with %s stop_on_gene can't be written : gene_len : %s" %
                          (exon.gene.name, exon.position, exon.stop, exon.gene.length))
                    count += 1
    print("total missed %s : %s" % (splicing_site, count))


def rnafold_launcher(rnafold_launcher, input_file, output_file):
    """
    Launch RNAfold.

    :param rnafold_launcher: (string) path where rnafold is installed
    :param input_file: (string) input fasta file
    :param output_file: (string) output rnafold.
    """
    subprocess.check_call("%sRNAfold --noPS %s > %s" % (rnafold_launcher, input_file, output_file), shell=True)


def extract_mfe(filename):
    """
    Extract the list of minimum free energy (MFE) from a file ``filename`` corresponding to \
    an output of RNAfold
    :param filename: (string) the output of RNAfold
    :return: (list of float) list of MFE values
    """
    list_mfe = []
    with open(filename, "r") as outfile:
        for line in outfile:
            if "." in line:
                line = line.split(" ")[-1]
                line = line.replace("(", "").replace(")", "")
                list_mfe.append(float(line))
    return list_mfe


def mfe_calculator(exon_list):
    """
    Launch svm_bp finder for every exon in the exon list
    :param exon_list: (list of ExonClass object) list of exons
    :return: (2 lists of floats) lists of float
    """
    output = os.path.realpath(os.path.dirname(__file__)) + "/control_dictionaries/"
    rnafold = os.path.realpath(os.path.dirname(__file__)).replace("src/minimum_free_energy", "data/RNAfold/")
    if not os.path.isdir(output):
        os.mkdir(output)
    file_3ss = "%sinput_rnafold_3ss.fa" % output
    out_3ss = "%sout_rnafold_3ss.fa" % output
    file_5ss = "%sinput_rnafold_5ss.fa" % output
    out_5ss = "%sout_rnafold_5ss.fa" % output
    fasta_writer(exon_list, file_3ss, "3ss")
    fasta_writer(exon_list, file_5ss, "5ss")
    rnafold_launcher(rnafold, file_3ss, out_3ss)
    rnafold_launcher(rnafold, file_5ss, out_5ss)
    mfe_3ss = extract_mfe(out_3ss)
    mfe_5ss = extract_mfe(out_5ss)
    return mfe_3ss, mfe_5ss
