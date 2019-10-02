#!/usr/bin/env python3

# -*- coding: UTF-8 -*-

"""
Description:
    The goal of this script is to test the intersection between experimental \
    and predicted branch point bed files.
"""


import subprocess
import os


def compute_intersection(experimental_bed, predicted_bed, output):
    """
    Compute the intersection between experimental and predicted bed of \
    branch points.

    :param experimental_bed: (str) a bed file containing experimental branch \
    points.
    :param predicted_bed: (str) a bed file containing predicted branch points.
    :param output: (str) folder where the intersection file will be created.
    :return: (str) the name of the created file
    """
    subprocess.check_call("intersectBed -a %s -b %s -s -u > "
                          "%s/predicted_intersect.bed" %
                          (predicted_bed, experimental_bed, output),
                          shell=True, stderr=subprocess.STDOUT)
    return "%s/predicted_intersect.bed" % output


def read_line(my_file):
    """
    Return the number of line contained in ``my_file``.

    :param my_file: (str) a file
    :return: (int) the number of lines in ``my_file``
    """
    return len(open(my_file, "r").readlines())


def get_only_regulated_exon(predicted_bed, output):
    """
    Select only the lines with GC and AT exons
    :param predicted_bed: (str) a bed containing the predicted branch points.
    :param output: (str) folder where the result will be created
    :return: (str) a bed file containing only AT or GC branch points.
    """
    new_pred = "%s/predicted_AT-GC_bp.bed.bed" % output
    ouf = open(new_pred, "w")
    with open(predicted_bed, "r") as infile:
        for line in infile:
            if "GC-exons" in line or "AT-exons" in line:
                ouf.write(line)
    ouf.close()
    return new_pred


def write_intersection(predicted_bed, experimental_bed, output, open_mode="w"):
    """
    Write the intersection file between ``predicted_bed`` and \
    ``experimental_bed``.

    :param predicted_bed: (str) a bed containing predicted branch points
    :param experimental_bed: (str) a bed containing experimental branch points.
    :param output: (str) folder where the result will be created
    :param open_mode: (str) the opening mode
    """
    intersect_bed = compute_intersection(predicted_bed, experimental_bed,
                                         output)
    with open("%s/intersection_predicted_exp_result.txt" % output, open_mode) \
            as ouf:
        nb_pred = read_line(predicted_bed)
        intersection_nb = read_line(intersect_bed)
        ouf.write("Predicted bed : %s : %s bp\n"
                  % (predicted_bed, nb_pred))
        ouf.write("experimental bed : %s : %s bp\n"
                  % (experimental_bed, read_line(experimental_bed)))
        ouf.write("intersection : %s : %s bp -> %s percent common bp with "
                  "prediction\n\n\n" % (intersect_bed, intersection_nb,
                                    round(intersection_nb / nb_pred * 100, 2)))



def main():
    """
    Test the size of the intersection between experimental \
    and predicted branch point bed files
    """
    base = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    output = base + "/result/experimental_branch_point"
    predicted_bp_bed = output + "/predicted_branch_points.bed"
    at_gc_pred = get_only_regulated_exon(predicted_bp_bed, output)
    exp_bp_bed = output + "/merged_branch_point_all_intron_filter.bed"
    write_intersection(predicted_bp_bed, exp_bp_bed, output, "w")
    write_intersection(at_gc_pred, exp_bp_bed, output, "a")


if __name__ == "__main__":
    main()
