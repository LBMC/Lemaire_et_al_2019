#!/usr/bin/env python3

# -*- coding: utf-8 -*-


"""
Description:
    The goal of this script is to create a table combining every info of interesert
"""
import group_factor
import union_dataset_function
import pandas as pd
import sqlite3
import control_exon_adapter
import os
import argparse
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3, "S": 4, "W": 5, "R": 6, "Y": 7}


def get_ase_events(cnx, id_project, regulation):
    """
    Get every exon up or down regulated in a particular project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :param regulation: (string)) up or down
    :return: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    """
    if regulation == "up":
        regulation = ">= 0.1"
    else:
        regulation = "<= -0.1"
    cursor = cnx.cursor()
    query = """SELECT gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue_glm_cor <= 0.05""" % (id_project, regulation)
    cursor.execute(query)
    res = cursor.fetchall()
    if len(res) == 0:
            query = """SELECT gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = %s
               AND delta_psi %s
               AND pvalue <= 0.05""" % (id_project, regulation)
            cursor.execute(query)
            res = cursor.fetchall()
    return res

def get_list_of_value_iupac_dnt(cnx, exon_list, target_column, nt_dnt):
    """
    Get the individual values of nt ``nt`` in ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :param nt_dnt: (string) a nucleotide or di_nucleotide
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            if r is not None:
                res.append(float(r.split(";")[nt_dic[nt_dnt]]))
            else:
                res.append(None)
    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                if r is not None:
                    res.append(float(r.split(";")[nt_dic[nt_dnt]]))
                else:
                    res.append(None)
                redundancy_gene_dic[exon[0]] = 1
    return res

def get_list_of_value(cnx, exon_list, target_column):
    """
    Get the individual values for ``target_column`` of every exon in ``exon_list``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param exon_list: (list of tuple of 2 int) each sublist corresponds to an exon (gene_id + exon_position on gene)
    :param target_column: (string) the column for which we want to get information on exons.
    :return: (list of float) values of ``target_column`` for the exons in  ``exon_list``.
    """
    cursor = cnx.cursor()
    res = []
    if target_column not in ["gene_size", "nb_intron_gene", "median_intron_size", "iupac_gene", "dnt_gene"]:
        for exon in exon_list:
            query = """SELECT %s
                       FROM sed
                       where gene_id = %s
                       AND exon_pos = %s """ % (target_column, exon[0], exon[1])
            cursor.execute(query)
            r = cursor.fetchone()[0]
            res.append(r)

    else:
        redundancy_gene_dic = {}
        for exon in exon_list:
            if exon[0] not in redundancy_gene_dic.keys():
                query = """SELECT %s
                           FROM sed
                           where gene_id = %s
                           AND exon_pos = %s """ % (target_column, exon[0], exon[1])
                cursor.execute(query)
                r = cursor.fetchone()[0]
                res.append(r)
                redundancy_gene_dic[exon[0]] = 1
    return res


def get_values_for_many_projects_iupac_dnt(cnx, id_projects_sf_names, target_columns, regulation, ctrl_full,
                                           exon_type):
    """
    Return the frequency of the nucleotide ``nt`` of ``target_column`` for each ``regulation`` \
    exons for projects in ``id_projects``.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_projects_sf_names: (list of str or  int) list project id if union is none. List of sf_name \
    else
    :param target_column: (string) the column for which we want to get information on exons.
    :param regulation: (string) up or down
    :param nt_dnt: (string) a nucleotide or a di-nucleotide
    :param union: (None or string) None if we want to work project by project, anything else to work \
    with exons regulation by a particular splicing factor.
    :return: (list of list of float) each sublist of float corresponds to the values of ``target_column`` \
    for every regulated exon in a given project.
    """
    print(target_columns)
    results = {target_column: [] for target_column in target_columns}
    results["project"] = []
    for sf_name in id_projects_sf_names:
        exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
        for target_column in target_columns:
            if "_nt_" in target_column:
                my_nt, target_column_name = target_column.split("_nt_")
                target_column_name = "iupac_" + target_column_name
                results[target_column] += get_list_of_value_iupac_dnt(cnx, exon_list, target_column_name, my_nt)
            else:
                results[target_column] += get_list_of_value(cnx, exon_list, target_column)
        results["project"] += [sf_name] * len(exon_list)
        # print("SF : %s" % sf_name)
        # for target_column in target_columns:
        #     print("%s len : %s"  % (target_column, len(results[target_column])))
        # print("project len : %s" % len(results["project"]))
    for target_column in target_columns:
        print(target_column)
        if "_nt_"in target_column:
            my_nt, target_column_name = target_column.split("_nt_")
            target_column_name = "iupac_" + target_column_name
            print("test len %s : %s" % (target_column_name, len(ctrl_full[target_column_name][my_nt])))
            results[target_column] += ctrl_full[target_column_name][my_nt]
            my_len = len(ctrl_full[target_column_name][my_nt])

        else:
            results[target_column] += ctrl_full[target_column]
            my_len = len(ctrl_full[target_column])
    results["project"] += [exon_type] * my_len
    # print("CCE")
    # print(ctrl_full.keys())
    # print(len(ctrl_full["iupac_upstream_intron"]["S"]))
    # print(type(ctrl_full["iupac_upstream_intron"]["S"]))
    for target_column in target_columns:
        print("%s len : %s" % (target_column, len(results[target_column])))
    print("project len : %s" % len(results["project"]))
    df = pd.DataFrame(results)
    return df



def main(target_columns, nt_list, output_folder, seddb, exon_type, size_bp_up_seq, regulation, name_tab):
    if seddb is None:
        seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    cnx = sqlite3.connect(seddb)
    ctrl_dic, ctrl_full = control_exon_adapter.control_handler(cnx, exon_type, size_bp_up_seq)
    sf_names = group_factor.get_wanted_sf_name(cnx)
    target_columns_new = [target_columns[i].replace("iupac", "%s_nt" % nt_list[i]) for i in range(len(target_columns))]
    df = get_values_for_many_projects_iupac_dnt(cnx, sf_names, target_columns_new, regulation, ctrl_full,
                                           exon_type)
    if output_folder:
        df.to_csv("%s/%s" % (output_folder, name_tab), sep="\t", index=False)
    return df


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Create a heatmap with  the wanted parameters

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("required argument")


    required_args.add_argument('--columns', dest='columns',
                               help="The columns (in sed database) you want to analyze, they must be coma separated",
                               required=True)
    required_args.add_argument('--nt_list', dest='nt_list',
                               help="The list of nucleotide you want to analyse, the length of the list must be equal to the length of columns argument. nt must be coma separated",
                               required=True)

    parser.add_argument("--output", dest="output", help="folder_were_ the result will be created", default=".")
    parser.add_argument("--exon_type", dest="exon_type", help="The type of control exons", default="CCE")
    parser.add_argument("--size", dest="size", help="The size of upstream exon selected to find some features such as nb_branch_point", default=100)
    parser.add_argument("--regulation", dest="regulation",
                        help="The regulation wanted (up or down)",
                        default="down")
    parser.add_argument("--name", dest="name",
                        help="The name of the output table",
                        default="result_tab.txt")
    required_args.add_argument("--seddb", dest="seddb", help="file corresponding to sed database", required=True)

    args = parser.parse_args()  # parsing arguments


    args.columns = args.columns.split(",")
    args.nt_list = args.nt_list.split(",")
    if len(args.columns) != len(args.nt_list):
        parser.error("ERROR : the lenght of the arguments columns and nt_list differt !")
    if not os.path.isdir(args.output):
        parser.error("ERROR;: the output folder does not exist")
    if args.exon_type not in ["CCE", "ACE", "ALL"]:
        parser.error("ERROR: argument exon_type must equals to CCE, ACE or ALL")
    try:
        args.size = int(args.size)
        if args.size not in [25, 50 ,35, 100]:
            parser.error("ERROR : the size must be equals to 25, 50, 35 or 100")
    except ValueError:
        parser.error("ERROR : the size must be an integer in 25, 50, 35 or 100")

    if args.regulation not in ["down", "up"]:
        parser.error("ERROR : wrong argument regulation : it must be up or down")

    main(args.columns, args.nt_list, args.output, args.seddb, args.exon_type, args.size, args.regulation, args.name)  # executing the program


if __name__ == "__main__":
    launcher()