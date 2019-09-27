#!/usr/bin/env python3

"""
Description:
    The goal of this script is to merge every branch point present \
    in a list of branch point bed file to remove redundancy
"""

import os
import sys
import subprocess as sub
import sqlite3


def read_bed(bed_file):
    """
    Get every line of a bed file.

    :param bed_file: (str) a bed file
    :return: (list of list of data) every line in the bed file
    """
    bed_list = []
    with open(bed_file, "r") as bedin:
        for line in bedin:
            if line[0] != "#":
                line = line.replace("\n", "")
                line = line.split("\t")
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_list.append(line)
    return bed_list


# def has_near(dic, dist_near, tcoordinates):
#     """
#     Find if the interval ``coordinates`` is near to another one.
#
#     :param dic: (dictionary) a dictionary of bed lines
#     :param dist_near: (int) the distance to detect the near branch point
#     :param tcoordinates: (list of 5 values) a genomic interval
#     :return: (boolean) True if the interval has near branch point false else
#     """
#     for bp in dic.keys():
#         tbp = dic[bp]
#         if tbp["coord"][0] == tcoordinates[0] and \
#             tbp["coord"][3] == tcoordinates[3]:
#             coordinates = tcoordinates.copy()
#             coordinates[1] -= dist_near
#             coordinates[2] += dist_near
#             if coordinates[1] <= tbp["coord"][1] < coordinates[2]:
#                 # print(tbp)
#                 return True
#     return False


def merge_bed_file(bed1, dic_bed):
    """

    :param bed1: (str) a bed file
    :param dic_bed: (dictionary) a dictionary of bed lines
    :return: (dictionary) dic_bed updated
    """
    bed_lines = read_bed(bed1)
    count_near = 0
    dist_near = 5
    tot = len(bed_lines)
    mcount = 0
    for bed_line in bed_lines:
        coordinates = [bed_line[0], bed_line[1], bed_line[2], bed_line[5]]
        key_name = ":".join(map(str, coordinates))
        if key_name not in dic_bed.keys():
            # res = has_near(dic_bed, dist_near, coordinates)
            # if res:
            #     # print(bed_line)
            #     # exit(1)
            #     count_near += 1
            dic_bed[key_name] = {"line": bed_line, "coord": coordinates}
        mcount += 1
        sys.stdout.write("%s / %s    --  %s       \r" % (mcount, tot, count_near))
    print(" --> found %s near exons for %s exon added in a dic of %s lines"
          % (count_near, len(bed_lines), len(dic_bed.keys())))
    return dic_bed


def write_bed(output, dic_bed, name):
    """
    Write the bed file.

    :param output: (str) folder where the result will be created
    :param dic_bed:  (dictionary) a dictionary of bed lines
    """
    filename = "%s/%s.bed" % (output, name)
    with open(filename, "w") as ouf:
        for bp in dic_bed.keys():
            ouf.write("\t".join(map(str, dic_bed[bp]["line"])) + "\n")
    return filename


def prepare_bed(folder, chrom_size_file, name):

    add_5 = "{0}/{1}_add_5.bed".format(folder, name)
    cmd = "bedtools slop -i {0}/{1}.bed -b 5 -g {2} > {3}".format(folder, name, chrom_size_file, add_5)
    sub.check_call(cmd, shell=True, stderr=sub.STDOUT)
    sorted_bed = add_5.replace(".bed", "_sorted.bed")
    cmd2 = "bedtools sort -i {0}  > {1}".format(add_5, sorted_bed)
    sub.check_call(cmd2, shell=True, stderr=sub.STDOUT)
    merged_bed = sorted_bed.replace(".bed", "_merged.bed")
    cmd3 = "bedtools merge -i {0} -s -d -1 > {1}".format(sorted_bed, merged_bed)
    sub.check_call(cmd3, shell=True, stderr=sub.STDOUT)


def get_intron_coordinates(cnx, output):
    """
    Write an intron proxi bed.

    :param cnx: (sqlite3 connectio object) connection to fasterdb lite
    :return:
    """
    dic = {1: "+", -1: "-"}
    count = 0
    cursor = cnx.cursor()
    query = """
    SELECT t2.id_gene, t2.pos_on_gene, t1.chromosome, 
           t2.start_on_chromosome - 1, t2.end_on_chromosome, 
           t1.strand
    FROM genes t1, introns t2
    WHERE t1.id = t2.id_gene
    """
    content_intron_bed = []
    cursor.execute(query)
    result = cursor.fetchall()
    for intron in result:
        intron = list(intron)
        if intron[3] < intron[4]:
            intron[5] = dic[intron[5]]
            if intron[5] == "+" and intron[4] - intron[3] > 100:
                intron[3] = intron[4] - 100
            if intron[5] == "-" and intron[4] - intron[3] > 100:
                intron[4] = intron[3] + 100
            my_row = ["chr" + intron[2], intron[3], intron[4]] + \
                     ["%s_%s" % (intron[0], intron[1])] + \
                     [".", intron[5]]
            content_intron_bed.append("\t".join(list(map(str, my_row))))
        else:
            count += 1
    print("Exon with a negative or nul size : %s" % count)
    filename = "%s/intron_proxi.bed" % output
    with open(filename, "w") as ouf:
        ouf.write("\n".join(content_intron_bed) + "\n")
    return filename


def filter_on_introns_proxi(merged_bed, intron_bed, output, final_name):
    """
    Create a bed of merged branch point only falling in ``intron_bed``.

    :param merged_bed: (str) a bed file containing branch point
    :param intron_bed: (str) a bed file containing the 100 last nt of introns
    :param output: (str) the output folder
    :param final_name: (str) the name of the final file
    """
    final_name = "%s/%s" % (output, final_name)
    cmd = "intersectBed -a %s -b %s -s -u > %s" % \
          (merged_bed, intron_bed, final_name)
    sub.check_call(cmd, shell=True, stderr=sub.STDOUT)


def filter_bed(my_bed, threshold, intron_bed, output):
    """
    :param my_bed: (str) a bed file
    :param threshold: (int) a threshold
    :param intron_bed: (str) a bed of the intronic region near the 3'ss
    :param output: (str) folder where results will be created
    :return: (str) the new name of the bed
    """
    if "mercer" in my_bed:
        new_bed = my_bed
    else:
        new_name = my_bed.replace(".bed", "_%scov_filter.bed" % threshold)
        new_bed = "%s/%s" % (output, os.path.basename(new_name))
        ouf = open(new_bed, "w")
        with open(my_bed, "r") as infile:
            for line in infile:
                line = line.replace("\n", "")
                if int(line.split("\t")[4]) >= threshold:
                  ouf.write(line + "\n")
        ouf.close()
    final_bed = new_bed.replace(".bed", "_intron_end_filter.bed")
    filter_on_introns_proxi(new_bed, intron_bed, output,
                            os.path.basename(final_bed))
    return final_bed


def main():
    """
    Create the merged bed of branch points.
    """
    threshold = 5
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__))))
    fasterdb = base_dir + "/data/fasterDB_lite.db"
    output = base_dir + "/result/experimental_branch_point"
    branch_point_folder = base_dir + "/data/experimental_branch_point"
    bp_bed_list = [branch_point_folder + "/" + mfile for mfile in
                   os.listdir(branch_point_folder) if "~" not in mfile
                   and mfile[0] != "."]
    if not os.path.isdir(output):
        os.mkdir(output)
    cnx = sqlite3.connect(fasterdb)
    intron_bed = get_intron_coordinates(cnx, output)
    dic_bed = {}
    list_filtered_bed = []
    for bed_file in bp_bed_list:
        print("Working on %s" % bed_file)
        dic_bed = merge_bed_file(bed_file, dic_bed)
        list_filtered_bed.append(filter_bed(bed_file, threshold, intron_bed,
                                            output))
    mname = "merged_branch_point_all"
    merged_bed = write_bed(output, dic_bed, mname)
    final_bed = "merged_branch_point_all_intron_filter.bed"
    filter_on_introns_proxi(merged_bed, intron_bed, output, final_bed)
    cnx.close()

    # merge taggart and pineda filtered
    dic_bed = {}
    for bed_file in list_filtered_bed:
        if "mercer" not in bed_file:
            dic_bed = merge_bed_file(bed_file, dic_bed)
    mname = "pineda-taggart_filtered_merged"
    write_bed(output, dic_bed, mname)


if __name__ == "__main__":
    main()

