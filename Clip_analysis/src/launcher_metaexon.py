#!/usr/bin/env python3

import lazyparser as lp
import os
from metaexon_coverage_ctrl import find_clip_bed_files, coverage_wrapper


def find_wanted_exon_list(folder):
    """
    Get the up and down regulated exons files by ``sfname`` in ``folder``.

    :param folder: (str) a directory containing exon file
    :return: (list of 2 values) the list of wanted exon list.
    """
    mfiles = sorted(os.listdir(folder))
    mfiles = ["%s/%s" % (folder, mfile) for mfile in mfiles]
    names = [os.path.basename(mfile).split(".")[0] for mfile in mfiles]
    if len(mfiles) != 2:
        raise IndexError("Error the len of mfiles is %s while it should be 2" %
                         len(mfiles))
    return mfiles, names


@lp.parse(input_folder="dir", folder_exon="dir", exon_bed="file",
          chrom_size_file="file", output="dir", metagene_script="file")
def main(input_folder, folder_exon, exon_bed, chrom_size_file, output,
         metagene_script):
    """

    :param input_folder: (str) a folder containing clip bed files files
    :param folder_exon: (str) a folder containing exon list
    :param exon_bed: (str) a bed file corresponding to fasterdb exons
    :param chrom_size_file: (str) a file indicating the length of every \
    human chromosome
    :param output: (str) path where the result will be created
    :param metagene_script: (str) path to the metagene script.
    """
    clip_beds = find_clip_bed_files(input_folder)
    exon_files, exonf_names = find_wanted_exon_list(folder_exon)
    color = "00aa00,5555FF"
    if exon_files[0] is not None:
        for clip_bed in clip_beds:
            sfname = os.path.basename(clip_bed).split("_")[0].upper()
            print("Working on %s : files => %s" % (sfname, clip_bed))
            for pos in ["n", "n-1", "n+1"]:
                coverage_wrapper(clip_bed, exon_bed, exon_files, exonf_names,
                                 chrom_size_file, output, metagene_script, "",
                                 color, "", pos)
    else:
        print("Error exon_files is not defined")
        exit(1)


if __name__ == "__main__":
    main()
