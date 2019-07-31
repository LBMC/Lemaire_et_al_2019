#!/usr/bin/env python3

"""
The goal of this script is to find good GC candidate exons to validate the \
model. This script will return every GC exons regulated by SNRPC, SNRNP70 and \
DDX5-D17 with the size of their flanking introns, the strength of \
the 5' splicing site and their minimum free energy
"""

import sqlite3
import subprocess
import lazyparser as lp
import sys
import pandas as pd
import os
from io import StringIO
from functools import reduce
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
import union_dataset_function as udf
fold = os.path.abspath(os.path.dirname(os.path.dirname(__file__))) + \
    "/minimum_free_energy/"
sys.path.insert(0, fold)
import function
import exon_class
fold_bp = os.path.abspath(os.path.dirname(os.path.dirname(__file__))) + \
    "/make_control_files_bp_ppt/"
sys.path.insert(0, fold_bp)
import function_bp
import exon_class_bp


def get_exon(exon_file):
    """
    Return the list of exon of interest.

    :param exon_file: (str) a file containing exons
    :return: (list of list of 2 int) the list of exons of interest
    """
    list_exon = []
    with open(exon_file, "r") as infile:
        for line in infile:
            list_exon.append(list(map(int,
                                      line.replace("\n", "").split("\t"))))
    return list_exon


def get_union_exon(exon_list1, exon_list2):
    """
    Return the exons in ``exon_list1`` and ``exon_list2``

    :param exon_list1: (list of list of 2 int) exons list 1
    :param exon_list2: (list of list of 2 int) exons list 2
    :return: (list of list of 2 int) the exon shared between ``exon_list1`` \
    and ``exon_list2``
    """
    return [exon for exon in exon_list1 if exon in exon_list2]


def get_exon_data(cnx, exon_list, ss):
    """
    Get the size of flanking introns + the strength of 5'SS for \
    every exons in ``exon_list``.

    :param cnx: (sqlite3 connection objects) connection to sed database
    :param exon_list: (list of list of 2 int) list of exons
    :param ss: (str) the splicing site of interest
    :return: (pandas dataframe)
    """
    if ss == "5'ss":
        force = "force_donor"
    else:
        force = "force_acceptor"
    res = []
    cursor = cnx.cursor()
    for exon in exon_list:
        query = "SELECT gene_symbol, gene_id, exon_pos, upstream_intron_size, " \
                "downstream_intron_size, %s " \
                "FROM sed " \
                "WHERE gene_id = %s and exon_pos = %s" % \
                (force, exon[0], exon[1])
        cursor.execute(query)
        tmp = cursor.fetchall()
        if len(tmp) > 1:
            raise IndexError("Only one exons should be found for %s_%s" %
                             (exon[0], exon[1]))
        res.append(list(tmp[0]))
    df = pd.DataFrame(res, columns=["gene_name", "gene_id", "pos",
                                    "upstream_intron_size",
                                    "downstream_intron_size",
                                    force])
    return df


def computing_mfe(cnx, df, output):
    """
    Add a column mfe 5'ss to the existing dataframe.

    :param cnx: (sqlite3 dataframe object) connection to fasterdb
    :param df: (pandas dataframe) table of exons
    :param output: (str) files were the mfe results will be created
    :return:  (pandas dataframe) table of exons with mfe data
    """
    exon_class_list = []
    exon_list = df[["gene_name", "gene_id", "pos"]].values
    for exon in exon_list:
        exon_class_list.append(exon_class.ExonClass(cnx, exon[0],
                                                    exon[1], exon[2]))
    mfe_3ss, mfe_5ss = function.mfe_calculator(exon_class_list,
                                               output, ps=True)
    df["mfe_5ss"] = mfe_5ss
    return df


def run_svs_bp_finder(input_file, gene_id, pos, sequence):
    """
    Run svm bp inder on ``input_file`` fasta file input.

    :param input_file: (str) a fasta file
    :param gene_id: (int) the id of the gene
    :param pos: (int) the pos of the exon
    :param sequence: (str) the intronic sequence
    :return: (pandas dataframe)
    """
    svm_prog = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    svm_prog += "data/SVM_BP_finder/svm_bpfinder.py"
    result = subprocess.check_output([svm_prog, "--input", input_file,
                                      "--species", "Hsap", "-l",
                                      "100"]).decode("ascii")
    result = StringIO(result)
    df = pd.read_csv(result, sep="\t")
    subprocess.check_call(['rm', input_file])
    ttt_motifs = sequence.count("TTT")
    df = df[(df["svm_scr"] > 0)]
    if df.empty:
        df["gene_id"] = []
        df["pos"] = []
        df["3'ss_intron_seq"] = []
        df["nb_TTT_motifs"] = []
        df.append({c: "NA" for c in df.columns}, ignore_index=True)

    else:
        df["bp_seq"] = df["bp_seq"].apply(lambda x: x.upper())
        df["gene_id"] = [int(gene_id)] * df.shape[0]
        df["pos"] = [int(pos)] * df.shape[0]
        df["3'ss_intron_seq"] = [sequence] * df.shape[0]
        df["nb_TTT_motifs"] = [ttt_motifs] * df.shape[0]
    df = df[["gene_id", "pos", "ss_dist", "bp_seq", "y_cont", "svm_scr",
             "3'ss_intron_seq", "nb_TTT_motifs"]]
    return df


def svm_bp_finder_launcher(cnx, exon_list, output):
    """
    Compute the number of good branch points of every exons in ``exon_list``.

    :param cnx: (sqlite3 connect object) connection to fasterDB database
    :param exon_list: ((list of list of 1 str 2 int) list of exons.
    :param output: (str) folder were the input will be created
    :return: (pandas dataframe) list of the branch point of interest
    """
    list_df = []
    for exon in exon_list:
        class_exon = exon_class_bp.ExonClass(cnx, exon[0], exon[1], exon[2])
        input_file = function_bp.fasta_writer(class_exon, output, 100)
        df = run_svs_bp_finder(input_file, exon[1], exon[2],
                               class_exon.upstream_intron.sequence_proxi)
        list_df.append(df)
    return pd.concat(list_df, ignore_index=True)



@lp.parse(exon_file="file", sed="file", fasterdb="file", output="dir",
          ss=["5'ss", "3'ss"])
def main(exon_file, name_table, list_sf, sed, fasterdb, output,
         ss="5'ss"):
    """
    Create a table showing for the exon commons in exon_files files \
    their surrounding introns length and their MFE at their 5'ss.

    :param exon_file: (str) a file containing gc/at exons
    :param name_table: (str) the name of the resulting table
    :param list_sf: (List(vtype=str)) list of sf name
    :param sed: (str) path to sed database
    :param fasterdb: (str) path to fasterdb database
    :param output: (str) file were the output will be created
    :param ss: (str) the splicing site of interest
    """
    sf_names = "_".join([name_table] + list_sf)
    exon_class.set_debug(1)
    exon_class_bp.set_debug(debug=1)
    cnx_sed = sqlite3.connect(sed)
    cnx_fasterdb = sqlite3.connect(fasterdb)
    exon_list = []
    print("Getting exon from file")
    exon_list.append(get_exon(exon_file))
    print("Getting regulated exons")
    for sf in list_sf:
        tmp = udf.get_every_events_4_a_sl(cnx_sed, sf, "down")
        tmp = [[int(v[0]), int(v[1])] for v in tmp]
        exon_list.append(tmp)
        print("\t%s : %s down-regulated exons" % (sf, len(tmp)))
    new_exon_list = reduce(get_union_exon, exon_list)
    print("Commons exons : %s" % len(new_exon_list))
    print("Getting commons exons data !")
    df = get_exon_data(cnx_sed, new_exon_list, ss)
    if ss == "5'ss":
        noutput = output + "/rnafold_" + sf_names + "_commons_down_exons/"
        print("Computing MFE")
        df = computing_mfe(cnx_fasterdb, df, noutput)
    else:
        # Code to compute number of good branch point
        print("Computing Good branch point")
        nexon_list = df[["gene_name", "gene_id", "pos"]].values
        df2 = svm_bp_finder_launcher(cnx_fasterdb, nexon_list, output)
        print(df2.head())
        print(df.head())
        df = pd.merge(df, df2, how="right", on=["gene_id", "pos"])
    print("Writing results !")
    df.to_csv("%s/%s_commons_down_exons.csv" % (output, sf_names), sep="\t", index=False)


if __name__ == "__main__":
    main()



