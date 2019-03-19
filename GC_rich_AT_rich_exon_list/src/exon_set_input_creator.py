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
import union_dataset_function

# functions
def connexion(seddb):
    """
    Connexion to SED database.

    :param seddb: ((string) path to sed database
    :return:  (sqlite3 connection object) allow connexion to sed database
    """
    return sqlite3.connect(seddb)


def get_splicing_factor_name():
    """
    Get the name of every splicing factor in splicing lore.

    :return: (list of string) list of splicing factor name
    """
    sf_name = ["DDX5_DDX17", "SNRNP70", "SNRPC", "U2AF1", "U2AF2", "SF1", "SF3A3", "SF3B4"]
    return sf_name


def input_writer_union(washed_exon_list, sf_name, output):
    """
    Write the input.

    :param washed_exon_list: (list of tuple of 2 int ans 1 str) each sublist corresponds to an exon \
    (gene_id + exon_position on gene + exon_regulation). Every exon regulated by a splicing factor in different \
    projects without redundancy.
    :param sf_name: (string) the name of the splicing factor regulating the exons in ``washed_exon_list``
    :param output: (string) path were the result will be created
    """
    with open("%sinput_%s-union.txt" % (output, sf_name.strip()), "w") as out_file:
        for exon in washed_exon_list:
            out_file.write(exon[0] + "\t" + exon[1] + "\n")


def main():
    """
    Execute the entire program
    """
    regulation = "down"
    output = os.path.realpath(os.path.dirname(__file__)).replace("src" , "result/inputs_union/")
    if not os.path.isdir(output):
        os.mkdir(output)
    seddb = output.replace( "result/inputs_union/", "data/sed.db")
    cnx = connexion(seddb)
    sf_names = get_splicing_factor_name()
    for sf_name in sf_names:
        print("Creating input for %s" % sf_name)
        exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
        input_writer_union(exon_list, sf_name, output)


if "__main__" == __name__:
    main()
