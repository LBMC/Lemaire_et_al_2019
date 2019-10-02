#!/usr/bin/env python3

# coding: utf-8

"""
Description:

    This script aims to create for each splicing factor an input (for ESA program) \
    corresponding to every exons regulated (always in the same way) by this splicing factor \
    in at least one cell line. The exons set in those input are called union-exons-set

Example:

    Let's say we have some sets of exons regulated by SRSF1 in two cell lines as follow:

    **SF:** SRSF1 - **Cell line:** Hela

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | SNRPC      |    3      |     UP        |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | MBNL1      |    9      |     DOWN      |
    +------------+-----------+---------------+


        **SF:** SRSF1 - **Cell line:** 293T

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | SNRPC      |    3      |     DOWN      |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | FUS        |    19     |     UP        |
    +------------+-----------+---------------+

    Then we will have the following exon set in result :

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | MBNL1      |    9      |     DOWN      |
    +------------+-----------+---------------+
    | FUS        |    19     |     UP        |
    +------------+-----------+---------------+

    .. note::

            * data about ``SNRPC_3`` is lost because it show different regulation in the cell lines of interest. \
            * ``hnRNPK_4`` exons is only displayed once (we do not repeat exon information). \
            * ``MBNL1_9`` and ``FUS_19`` are both displayed because they are regulated by SRSF1 in at least 1 cell line.

The input will of course respect ESA input format. See \
`ESA gitlab page <https://gitlab.biologie.ens-lyon.fr/Auboeuf/uniform_exons_features/Exon_Set_Analyzer>`_ for \
more information.
"""


# import
import sqlite3
import os
import sys


bad_id_projects = [139, 13, 164]
bad_factor = ["MBNL1_2", "PTBP1_2", "RBM10", "ESRP2",
              "RBM47", "RBM17", "TIA1", "SRSF10"]

# functions
def connexion(seddb):
    """
    Connexion to SED database.

    :param seddb: ((string) path to sed database
    :return:  (sqlite3 connection object) allow connexion to sed database
    """
    return sqlite3.connect(seddb)


def get_interest_project(cnx):
    """
    Get the id of every project defined in sed database (from splicing lore).

    :param cnx: (sqlite3 connection object) connexion to sed database
    :return: (list of int) list of id_project
    """
    cursor = cnx.cursor()
    query = "SELECT id, project_name FROM rnaseq_projects"
    cursor.execute(query)
    res = cursor.fetchall()
    idp = [val[0] for val in res]
    name = [val[1] for val in res]
    return idp, name


def get_splicing_factor_name(cnx):
    """
    Get the name of every splicing factor in splicing lore.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :return: (list of string) list of splicing factor name
    """
    cursor = cnx.cursor()
    query = "SELECT DISTINCT sf_name FROM rnaseq_projects"
    cursor.execute(query)
    res = cursor.fetchall()
    sf_name = [val[0] for val in res if val[0] not in bad_factor]
    return sf_name


def get_projects_links_to_a_splicing_factor(cnx, sf_name):
    """
    Get the id of every projects corresponding to a particular splicing factor.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_name: (string) splicing factor name
    :return: (list of int) list of id_project
    """
    cursor = cnx.cursor()
    query = "SELECT id FROM rnaseq_projects WHERE sf_name = ?"
    cursor.execute(query, (sf_name,))
    res = cursor.fetchall()
    idp = [val[0] for val in res if val[0] not in bad_id_projects]
    if sf_name == "SRSF9":
        print("------------------------SRSF9-----------------------")
        print(idp)
    return idp


def get_ase_events(cnx, id_project):
    """
    Get every exon regulated in a particular project.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param id_project: (int) a project id
    :return: (list of tuple of one str 2 int) each sublist corresponds to an exon (\
    exon_regulation + gene_id + exon_position on gene)
    """
    cursor = cnx.cursor()
    print("request1")
    query = """SELECT delta_psi, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = ?
               AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
               AND pvalue_glm_cor <= 0.05"""
    cursor.execute(query, (id_project,))
    print("handling request")
    res = cursor.fetchall()
    if len(res) == 0:
            print("request2")
            query = """SELECT delta_psi, gene_id, exon_skipped
               FROM ase_event
               WHERE id_project = ?
               AND (delta_psi >= 0.1 OR delta_psi <= -0.1)
               AND pvalue <= 0.05"""
            cursor.execute(query, (id_project,))
            res = cursor.fetchall()
    nres = []
    for exon in res:
        nexon = list(exon[1:3])
        if exon[0] < 0:
            nexon = ["down"] + nexon
        else:
            nexon = ["up"] + nexon
        nres.append(nexon)
    return nres


def get_every_events_4_a_sl(cnx, sf_name):
    """
    Get every splicing events for a give splicing factor.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_name: (string) the name of a splicing factor
    :return: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects.
    """
    exons_list = []
    print("getting id projects...")
    id_projects = get_projects_links_to_a_splicing_factor(cnx, sf_name)
    for id_project in id_projects:
        print("getting exon of projects %s" % id_project)
        ase_event = get_ase_events(cnx, id_project)
        exons_list += ase_event
    return exons_list


def washing_events(exon_list):
    """
    Remove redundant exons or remove exons showing different regulation.

    :param exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects.
    :return: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon (gene_id + exon_position on gene + \
    exon_regulation). Every exon regulated by a splicing factor in different projects without redundancy.
    """
    replace_dic = {"up":"down", "down": "up"}
    dic = {}
    prefix_list = []
    for exon in exon_list:
        exon_prefix = "%s_%s" % (exon[1], exon[2])
        exon_name = "%s_%s" % (exon[0], exon_prefix)
        if exon_name not in dic:
            if exon_prefix not in prefix_list:
                dic[exon_name] = 1
                prefix_list.append(exon_prefix)
            else:
                reverse_name = exon_name.replace(exon[0], replace_dic[exon[0]])
                if reverse_name in dic:
                    del(dic[reverse_name])
                # Else : the exon was deleted before because of a different regulation
        else:
            dic[exon_name] += 1
    # creation of the new list of exons
    new_exon_list = []
    for key in dic:
        my_exon = key.split("_")
        my_exon = ["_".join(my_exon[1:3])] + my_exon
        new_exon_list.append(my_exon)
    return new_exon_list


def input_writer_union(washed_exon_list, sf_name, output, regulation):
    """
    Write the input.

    :param washed_exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon \
    (gene_id + exon_position on gene + exon_regulation). Every exon regulated by a splicing factor in different \
    projects without redundancy.
    :param sf_name: (string) the name of the splicing factor regulating the exons in ``washed_exon_list``
    :param output: (string) path were the result will be created
    """
    with open("%s/input_%s-union_%s.txt" % (output, sf_name.strip(), regulation), "w") as out_file:
        for exon in washed_exon_list:
            if regulation in exon:
                out_file.write("_".join(exon[0].split("_")) + "\n")


def input_writer(exon_list, project_name, output, regulation):
    """
    Write the input.

    :param exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon \
    (gene_id + exon_position on gene + exon_regulation).
    :param project_name: (string) of the project regulating the exon in exon_list
    :param output: (string) path were the result will be created
    """
    with open("%sinput_%s_%s.txt" % (output, project_name.strip(), regulation), "w") as out_file:
        for exon in exon_list:
            if regulation in exon:
                exon = list(map(str, exon))
                exon = [exon[1] + "_" + exon[2]] + exon
                out_file.write("_".join(exon[0].split("_")) + "\n")


def main():
    """
    Execute the entire program
    """
    regulation = "down"
    output = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    output += "/data/input"
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = output.replace( "input", "sed.db")
    cnx = connexion(seddb)
    sf_names = get_splicing_factor_name(cnx)
    for sf_name in sf_names:
        print("Creating input for %s" % sf_name)
        exon_list = get_every_events_4_a_sl(cnx, sf_name)
        washed_exon_list = washing_events(exon_list)
        input_writer_union(washed_exon_list, sf_name, output, regulation)


if "__main__" == __name__:
    main()