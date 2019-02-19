#!/usr/bin/python3

"""
Description:

    This file contains every group of factor we want to study.
"""

# list of splicing factors
splicing_factors = ("PCBP2", "HNRNPA1", "HNRNPU", "QKI", "PTBP1",
                "TRA2A_B", "KHSRP", "MBNL1",
                "HNRNPL", "HNRNPK", "SRSF7", "HNRNPA2B1", "SFPQ",
                "RBM15", "HNRNPM", "FUS",
                "DAZAP1", "RBM39", "SRSF9", "RBM25", "RBM22", "HNRNPF", "SRSF5",
                "PCBP1", "RBFOX2", "HNRNPH1", "RBMX", "SRSF6", "MBNL2", "SRSF1",
                "SRSF2", "HNRNPC", "SRSF3")

# list of factors that compose U1 snrnp
u1_factors = ("SNRNP70", "SNRPC", "DDX5_DDX17")
# list of factor that compose U2 snrnp
u2_factors = ("SF1", "SF3A3", "SF3B1", "SF3B4", "U2AF1", "U2AF2")

other_spliceosome = ("PRPF8", "SNRNP200")

# The id project we have to delete
bad_id_projects = [139, 13, 164]
# 139 : PTBP1_ENCSR527IVX_K562, 13:  SF3B1_GSE65644_Hela, 164 :SRSF9_ENCSR113HRG_K562


# function
def get_projects_links_to_a_splicing_factor_list(cnx, sf_list):
    """
    Get the id of every projects corresponding to a particular splicing factor.

    :param cnx: (sqlite3 connection object) connexion to sed database
    :param sf_list: (list of strings) list of splicing factor name
    :return: (2 lists):
        * (list of int) the list of project we want to analye
        * (list of string) the list of projects name we want to  analyse
    """
    id_projects = []
    name_projects = []
    cursor = cnx.cursor()
    query = "SELECT id, sf_name, db_id_project, cl_name FROM rnaseq_projects WHERE sf_name = ?"
    for sf_name in sf_list:
        cursor.execute(query, (sf_name,))
        res = cursor.fetchall()
        for val in res:
            if val[0] not in bad_id_projects:
                id_projects.append(val[0])
                name_projects.append("%s_%s_%s" % (val[1], val[2], val[3]))
    return id_projects, name_projects


def get_id_and_name_project_wanted(cnx, sf_type):
    """
    Get the list of project of interest and the name of the project wanted
    :param cnx: (sqlite3 connect object) connection to sed database
    :param sf_type: (string) the type of sf we want to analyze
    :return: (2 lists):
        * (list of int) the list of project we want to analye
        * (list of string) the list of projects name we want to  analyse
    """
    good_sf = splicing_factors
    id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)

    return id_projects, name_projects


def get_wanted_sf_name(sf_type):
    """
    Return the list of splicing factors of interest.

    :param sf_type: (string) the type of sf wanted
    :return: (list of string) the list of sf of interest
    """

    name_projects = splicing_factors

    return name_projects