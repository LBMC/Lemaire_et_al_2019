#!/usr/bin/python3

"""
Description:

    This file contains every group of factor we want to study.
"""

# list of factors regulating at rich exons
at_rich_down = ("PCBP2", "HNRNPA1", "HNRNPU", "QKI", "PTBP1",
                "TRA2A_B", "KHSRP", "MBNL1",
                "HNRNPL", "HNRNPK", "SRSF7", "HNRNPA2B1", "SFPQ",
                "RBM15", "HNRNPM", "FUS",
                "DAZAP1", "RBM39")

# the list of SF tregulating gc rich genes
gc_rich_down = ("SRSF9", "RBM25", "RBM22", "HNRNPF", "SRSF5",
                "PCBP1", "RBFOX2", "HNRNPH1", "RBMX", "SRSF6", "MBNL2", "SRSF1")

other = ("SRSF2", "HNRNPC", "SRSF3")
other_spliceosome = ("U2AF1", "SF3B1")
# list of factors that compose U1 snrnp
u1_factors = ("SNRNP70", "SNRPC", "DDX5_DDX17")
# list of factor that compose U2 snrnp
u2_factors = ("SF1", "SF3A3", "SF3B4", "U2AF2")
chromatin_factors = ("DNMT3A", "EZH2", "KMT2A", "KMT2D", "MBD2", "MBD3", "SETD2", "SUV39H1",
                     "SUV39H2", "TDG", "TET2", "EED")

# -- additional list of exons --
ct_rich_down = ("SRSF9", "RBMX", "PCBP1", "HNRNPF", "SRSF5", "RBFOX2",
                "RBM25", "MBNL2", "RBM22", "HNRNPH1", "SRSF6", "SRSF2",
                "SRSF7", "DAZAP1", "HNRNPC", "HNRNPA1", "PCBP2", "QKI",
                "FUS", "HNRNPL", "HNRNPU")
ga_rich_down = ("SRSF1", "SFPQ", "MBNL1", "HNRNPM", "PTBP1", "KHSRP",
                "HNRNPA2B1", "RBM15", "RBM39", "HNRNPK", "TRA2A_B")


color_dic = {"GC_pure": "#5555FF", "AT_pure": "#00aa00", "GC_all": "#0394d9",
             "AT_all": "#03d994", "CCE": "red", "ALL": "red", "ACE": "red",
             "DDX5_DDX17": "#AA00FF", "DDX5_17": "#AA00FF", "SNRNP70": "#5555D0", "SNRPC": "cyan",
             "U1-factors": "#7777F0", "U2-factors": "#33EE33",
             "U2AF1": "#006400", "U2AF2": "olive", "SF1": "#55FF55", "SF3A3": "#D8EF48", "SF3B4": "#8FBC8F"}

color_dic_bright = {"GC_pure": "#AAAAFF", "AT_pure": "#99FF99", "GC_all": "#53D4FF",
                    "AT_all": "#50FFDF", "CCE": "#FF9999", "ALL": "#FF9999", "ACE": "#FF9999",
                    "DDX5_DDX17": "#DD99FF", "DDX5_17": "#DD99FF", "SNRNP70": "#9999F0", "SNRPC": "#77FFFF",
                    "U1-factors": "#9999FF", "U2-factors": "#66FF66",
                    "U2AF1": "#55AA55", "U2AF2": "#DDDD55", "SF1": "#AAFFAA", "SF3A3": "#E8FF88", "SF3B4": "#AFECAF"}

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
    if sf_type is None:
        good_sf = at_rich_down + gc_rich_down + other
        id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)
    elif sf_type == "GC_rich" or sf_type == "gc_rich_down":
        good_sf = gc_rich_down
        id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)
    elif sf_type == "AT_rich" or sf_type == "at_rich_down":
        good_sf = at_rich_down
        id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)
    elif sf_type == "CF" or sf_type == "chromatin_factors":
        good_sf = chromatin_factors
        id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)
    else:
        good_sf = u1_factors + u2_factors
        id_projects, name_projects = get_projects_links_to_a_splicing_factor_list(cnx, good_sf)
    return id_projects, name_projects


def get_wanted_sf_name(sf_type):
    """
    Return the list of splicing factors of interest.

    :param sf_type: (string) the type of sf wanted
    :return: (list of string) the list of sf of interest
    """
    if sf_type is None:
        name_projects = at_rich_down + gc_rich_down + other
    elif sf_type == "GC_rich" or sf_type == "gc_rich_down":
        name_projects = gc_rich_down
    elif sf_type == "AT_rich" or sf_type == "at_rich_down":
        name_projects = at_rich_down
    elif sf_type == "CF" or sf_type == "chromatin_factors":
        name_projects = chromatin_factors
    elif sf_type == "GA_rich" or sf_type == "ga_rich_down":
        name_projects = ga_rich_down
    elif sf_type == "CT_rich" or sf_type == "ct_rich_down":
        name_projects = ct_rich_down
    elif sf_type == "other":
        name_projects = other
    elif sf_type == "all":
        name_projects = at_rich_down + gc_rich_down + other + u1_factors + u2_factors + other_spliceosome
    else:
        name_projects = u1_factors + u2_factors
    return name_projects
