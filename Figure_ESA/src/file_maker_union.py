#!/usr/bin/env python3

# -*- coding : utf-8 -*-

import os
import sys
import union_dataset_function
import figure_producer
import numpy as np
import sqlite3

def main():
    regulation="down"
    target_column = sys.argv[1]
    seddb = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/data/sed.db"
    cnx = sqlite3.connect(seddb)
    output = "/".join(os.path.realpath(__file__).split("/")[:-2]) + "/result/%s_union_dataset.txt" % target_column
    # If the output directory does not exist, then we create it !
    name_projects = union_dataset_function.get_splicing_factor_name(cnx)
    outfile = open(output, "w")
    outfile.write("sf_regulation\tmedian_%s\n" % target_column)
    for sf_name in name_projects:
        exon_list = union_dataset_function.get_events_4_a_sl(cnx, sf_name, regulation)
        values = np.array(figure_producer.get_list_of_value(cnx, exon_list, target_column))
        median_obs = np.median(values[~np.isnan(values)])
        outfile.write("%s_%s\t%s\n" % (sf_name, regulation, median_obs))
    outfile.close()

if __name__ == "__main__":
    main()

