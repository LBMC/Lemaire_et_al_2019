#!/usr/bin/python3

"""
Description:

    This file contains every group of factor we want to study.
"""

# list of factors regulating at rich exons
at_rich_down = ("PCBP2", "HNRNPA1", "HNRNPU", "QKI", "SRSF3", "RBM17", "PTBP1",
                "TIA1", "HNRNPC",
                "TRA2A_B", "KHSRP", "MBNL1",
                "HNRNPL", "HNRNPK", "SRSF7", "HNRNPA2B1", "SFPQ",
                "RBM15", "HNRNPM", "FUS",
                "DAZAP1", "RBM39")

# the list of SF tregulating gc rich genes
gc_rich_down = ("SRSF9", "RBM25", "RBM22", "SRSF2", "HNRNPF", "SRSF5",
                "PCBP1", "RBFOX2", "HNRNPH1", "RBMX", "DDX5_DDX17", "SRSF1", "SRSF6", "MBNL2")

# list of factors that compose U1 snrnp
u1_factors = ("SNRNP70", "SNRPC")
# list of factor that compose U2 snrnp
u2_factors = ("SF1", "SF3A3", "SF3B1", "SF3B4", "U2AF1", "U2AF2")
chromatin_factors = ("DNMT3A", "EZH2", "KMT2A", "KMT2D", "MBD2", "MBD3", "SETD2", "SUV39H1",
                     "SUV39H2", "TDG", "TET2", "EED")