#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
This script will merge the six figures of every clip coverage figure produced
"""

import lazyparser as lp
import re
import os
import numpy as np
import subprocess


def check_prefixes(folder):
    """
    Get the prefixes of every files in ``folder``.

    :param folder: (str) a directory containing coverage files from clip \
    sequence.
    :return: (list of str) list of prefixes for files
    """
    figures = os.listdir(folder)
    figures = [fig for fig in figures if ".png" in fig]
    prefixes = [re.split(r"_n_|_n\+1_|_n-1_", fig)[0] for fig in figures]
    prefixes, count = np.unique(prefixes, return_counts=True)
    new_pref = [prefixes[i] for i in range(len(prefixes)) if count[i] == 6]
    if len(new_pref) == 0:
        raise IndexError("new_pref should not be empty !")
    return new_pref


def figure_merger(folder, prefix):
    """
    From a prefix merge every prossible figures.

    :param folder: (str) a folder where the coverage figure are located.
    :param prefix: (str) the prefix of the figure of interest
    """
    list_fig = ["%s/%s_%s.png" % (folder, prefix, cur) for cur in
                ["n-1_start", "n-1_stop", "n_start", "n_stop",
                 "n+1_start", "n+1_stop"]]
    cmd = "montage -geometry +1+1 -density 200 -compress jpeg -tile 6x1 " \
          "%s %s/%s_recap.pdf" % (" ".join(list_fig), folder, prefix)
    print(cmd)
    subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)


@lp.parse(folder="dir")
def main(folder):
    """
    Merge the coverage files produced by the metaexon_coverage.py script.

    :param folder: (str) a folder containing coverage figures
    """
    prefixes = check_prefixes(folder)
    for prefix in prefixes:
        print("Working on %s" % prefix)
        figure_merger(folder, prefix)


if __name__ == "__main__":
    main()
