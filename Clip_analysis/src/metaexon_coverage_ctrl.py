#!/usr/bin/env python3


"""
Description:
    The goal of this script is to create a figure showing the coverage of \
    a bed file in the form of a metaexon graphics.

    Requirements: betools and bedGraphToBigWig must be in your path.
"""

import subprocess
import os
import gzip
import lazyparser as lp


def chrom_name_adapter(bed_file, new_name):
    """
    Correct chromosome name.

    :param bed_file: (string) a bed obtain from a clip experiment
    :param new_name: (string) new name of the file
    :return: (string) the name of the
    """
    if "gz" in bed_file:
        bedin = gzip.open(bed_file, "rt")
    else:
        bedin = open(bed_file, "r")
    with open(new_name, "w") as bedout:
        for line in bedin:
            if line[0] == "#":
                bedout.write(line)
            else:
                line = line.split("\t")
                line[0] = line[0].replace("chr", "")
                if line[0] == "M":
                    line[0] += "T"
                bedout.write("\t".join(line))
    bedin.close()
    return new_name


def bed2bw(bed_file, refsize, outfolder):
    """
    Convert a bed file to a bigwig file.

    :param bed_file: (string) a bed file
    :param refsize: (string) a file containing the hg19 chromosome size.
    :param outfolder: (string) path where the bigwig file will be created.
    :return: (string) the name of the bigwig file
    """
    # adapting chromosome column.
    bed_chrom = outfolder + os.path.basename(bed_file)
    bed_chrom = bed_chrom.replace(".bed", "nochr.bed").replace(".gz", "")
    chrom_name_adapter(bed_file, bed_chrom)
    # sorting the bed
    sort_bed = bed_chrom.replace(".bed", ".sort.bed").replace(".gz", "")
    sort_bed = outfolder + os.path.basename(sort_bed)
    cmd_sort = "bedSort %s %s" % (bed_chrom, sort_bed)
    print(cmd_sort)
    subprocess.check_call(cmd_sort, shell=True, stderr=subprocess.STDOUT)
    # creation of a bed file
    bg_file = outfolder + os.path.basename(sort_bed).replace(".bed", ".bg")
    bw_file = bg_file.replace(".bg", ".bw")
    cmd_bg = "bedtools genomecov -bg -i %s -g %s > %s" % \
             (sort_bed, refsize, bg_file)
    print(cmd_bg)
    subprocess.check_call(cmd_bg, shell=True, stderr=subprocess.STDOUT)
    # transformation to bigwig
    cmd_bw = "bedGraphToBigWig %s %s %s" % (bg_file, refsize, bw_file)
    print(cmd_bw)
    subprocess.check_call(cmd_bw, shell=True, stderr=subprocess.STDOUT)
    return bw_file


def metagene_covergae_launcher(metagene_script, bw_file, cond, bed_annot,
                               annot_name, color, outname, outdir, off_set):
    """
    Launch the metagene script program.

    :param metagene_script: (string) path where the metagene script is located
    :param bw_file: (string) path to the bigwigfile
    :param cond: (string) a title
    :param bed_annot: (string) a bed file containing the chosen \
    annotation to work on for metagene creation
    :param annot_name: (string) the name of the annotation file
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    :param outdir:  (string) path where the result will be created
    :param off_set: (int) offset for 0 in the graphic

    """
    cmd = """Rscript %s -bw %s -cond %s -rep 1 -annot %s \
    -prefixes %s --color-pallette %s -sign %s \
    -out_dir %s -off_set %s -comp_pair %s"""
    cmd = cmd % (metagene_script, bw_file, cond, bed_annot, annot_name,
                 color, outname, outdir, off_set, cond)
    print(cmd)
    subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)


def extract_exon_list(filename):
    """
    Extract an exon list from a file named ``filename``.

    :param filename:  (string)  a file containing an exon list.
    :return: (list of 2 int) gene id and exon_pos
    """
    exon_list = []
    with open(filename, "r") as outfile:
        for line in outfile:
            line = line.replace("\n", "").split("\t")
            exon_list.append(line)
    return exon_list


def exon_reader(exon_files, pos):
    """
    Read an exon files
    :param exon_files: (str) an exon file containing exons id
    :param pos: (str) the position of the exon to consider (n-1, n or n+1)
    :return: (list of str) list of exon
    """
    with open(exon_files, 'r') as ef:
        exons = ef.read().splitlines()
    if pos != "n":
        new_exons = []
        for exon in exons:
            exon = exon.split("_")
            if pos == "n+1":
                exon[1] = str(int(exon[1]) + 1)
            else:
                exon[1] = str(int(exon[1]) - 1)
            new_exons.append("_".join(exon))
        return new_exons
    return exons


def handle_bed_creator(dic_bed, exon_lists, exonl_names, template_fodler,
                       chrom_size_file, pos):
    """
    Get the name of the annotation file.

    :param dic_bed: (dictionary of list) a bed dictionary with \
    the key corresponding to exons.
    :param exon_lists: (list of 2 list of string) list of exons
    :param exonl_names: (list of 2 string) list of the name of the exons lists
    :param template_fodler: (string) a folder containing template bed file.
    :param chrom_size_file: (string) a file containing the length of \
    hg19 chromosomes.
    :param pos: (str) the position of the exons to consider (n, n-1, n+1)
    :return: 4 elements:

        1. (string) the list of bed files (coma separated) corresponding to\
         intron-exon junctions containing the 3'SS
        2. (string) the list of bed files (coma separated) corresponding to\
         intron-exon junctions containing the 5'SS
        3. (string) the list of name (coma separated) for the bed files )\
         corresponding to intron-exon junctions \
        containing the 3'SS
        4. (string) the list of name (coma separated) for the bed files )\
         corresponding to intron-exon junctions\
         containing the 5'SS
    """
    final_template_start = []
    final_template_stop = []
    for i in range(len(exonl_names)):
        if not os.path.isfile("%s/%s_%s_exon_%s_add_i50_o200.bed" %
                              (template_fodler, exonl_names[i], pos, "start")):
            templates_bed = bed_creator(dic_bed, exon_lists[i],
                                        template_fodler, exonl_names[i],
                                        chrom_size_file, pos)
        else:
            start = "%s/%s_%s_exon_%s_add_i50_o200.bed" % \
                    (template_fodler, exonl_names[i], pos, "start")
            templates_bed = [start, start.replace("start", "stop")]
            print("recovering template files : %s" % start)
        final_template_start.append(templates_bed[0])
        final_template_stop.append(templates_bed[1])
    finale_names_start = []
    finale_names_stop = []
    for my_name in exonl_names:
        finale_names_start.append(my_name + "_start")
        finale_names_stop.append(my_name + "_stop")
    finale_names_start = ",".join(finale_names_start)
    finale_names_stop = ",".join(finale_names_stop)
    final_templates_start = ",".join(final_template_start)
    final_templates_stop = ",".join(final_template_stop)
    return final_templates_start, final_templates_stop, \
        finale_names_start, finale_names_stop


def load_chrom_size_dic(chrom_size_file):
    """
    Load in the form a a dictionary the chromosome size file.

    :param chrom_size_file: (float) chromosome size filz
    :return: (dictionary of int) dictionary that links each chromosome \
    to it's size
    """
    dic_size = {}
    with open(chrom_size_file, "r") as size_file:
        for line in size_file:
            line = line.replace("\n", "").split("\t")
            if "_" not in line[0]:
                dic_size[line[0]] = int(line[1])
    return dic_size


def bed_creator(dic_bed, exon_list, dest_folder, name_bed, chrom_size_file,
                pos):
    """
    Create a bed file containing all the exons regulated by ``sf_name`` \
    with ``regulation``.

    :param dic_bed: (dictionary of list) a bed dictionary with the \
    key corresponding to exons.
    :param exon_list: (list of string) 'gene id'_'exon pos'
    :param dest_folder: (string) path where the bed will be created
    :param name_bed: (string) name of the file
    :param chrom_size_file: (string) a file containing chromosome size
    :param pos: (str) the exons to consider
    :return: (string) the name of the bed file created
    """
    iexon = 50
    oexon = 200
    dic_size = load_chrom_size_dic(chrom_size_file)
    exon_info_left = []
    exon_info_right = []
    for exon in exon_list:
        if exon in dic_bed.keys():
            res = dic_bed[exon].copy()
            if res[3] == "-":
                res[3] = "-1"
            else:
                res[3] = "1"
            exon_left = None
            exon_right = None
            if res[1] - oexon >= 0 and res[1] + iexon <= dic_size[str(res[0])]:
                exon_left = [res[0], res[1] - oexon, res[1] + iexon] + \
                            [exon, ".", res[3]]
            if res[2] - iexon >= 0 and res[2] + oexon <= dic_size[str(res[0])]:
                exon_right = [res[0], res[2] - iexon, res[2] + oexon] + \
                             [exon, ".", res[3]]
            if res[3] == "-1":
                tmp = exon_left
                exon_left = exon_right
                exon_right = tmp
            if exon_left:
                exon_info_left.append("\t".join(list(map(str, exon_left))))
            if exon_right:
                exon_info_right.append("\t".join(list(map(str, exon_right))))
    contents = ["\n".join(exon_info_left), "\n".join(exon_info_right)]
    names = ["start", "stop"]
    final_names = []
    for i in range(len(contents)):
        filename = "%s/%s_%s_exon_%s_add_i50_o200.bed" % \
                   (dest_folder, name_bed, names[i], pos)
        final_names.append(filename)
        with open(filename, "w") as bedfile:
            bedfile.write(contents[i] + "\n")
    return final_names


def dic_bed_maker(bed):
    """
    Turn a bed into a dic bed.

    :param bed: (string) a bed file.
    :return: (dictionary of list) diciotnary of bed
    """
    dic = {}
    with open(bed, 'r') as bf:
        for line in bf:
            if line[0] != "#":
                line = line.replace("\n", "").split("\t")
                dic[line[3]] = [line[0], int(line[1]), int(line[2]), line[5]]
    return dic


def coverage_wrapper(clip_bed, exon_bed, exon_files, exonf_names, refsize,
                     output, metagene_script, cond, color, outname, pos):
    """

    :param clip_bed: (string) a bed file coming from a clip experiment
    :param exon_bed: (string) a bed corresponding to every fasterDB exons
    :param exon_files: (list of 2 string) a list of 2 files containing up and\
     down exons by the same splicing factor
    :param exonf_names: (list of 2 string) list of the name of the 2 exons \
    files
    :param refsize: (string) a file containing the hg19 chromosome size.
    :param output: (string) path where the result will be created.
    :param metagene_script: (string) path where the metagene script is located
    :param cond: (string) a title
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    :param pos: (str) the position of exons to consider (n, n-1, n+1)
    """
    if len(exon_files) != len(exonf_names):
        raise IndexError("Arguments exon_files and exonf_names \
                          don't have the same lenght !")
    template_fodler = os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)))
    template_fodler += "/result/.template"
    parent_template = os.path.dirname(template_fodler)
    if not os.path.isdir(os.path.dirname(template_fodler)):
        os.mkdir(parent_template)
    if not os.path.isdir(template_fodler):
        os.mkdir(template_fodler)
    dic_bed = dic_bed_maker(exon_bed)
    exon_lists = [exon_reader(mfile, pos) for mfile in exon_files]
    tplt_start, tplt_stop, n_start, n_stop = \
        handle_bed_creator(dic_bed, exon_lists, exonf_names,
                           template_fodler, refsize, pos)
    output_bw = output + "bigwig/"
    output_graph = output + "figure/"
    if not os.path.isdir(output_bw):
        os.mkdir(output_bw)
    if not os.path.isdir(output_graph):
        os.mkdir(output_graph)

    bw_file = bed2bw(clip_bed, refsize, output_bw)
    if not outname:
        outname = os.path.basename(clip_bed)
        outname = outname.replace(".gz", "").replace(".bed", "")
        outname += "_%s" % pos
    if not cond:
        cond = exonf_names[0].split("_")[0]
        cond += "_%s" % pos
    try:
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_start,
                                   n_start, color, outname + "_start",
                                   output_graph, 200)
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_stop,
                                   n_stop, color, outname + "_stop",
                                   output_graph, 50)
        print("\033[92m{0}SUCCESS ABOVE{0}\033[0m".format("-" * 30))
    except subprocess.CalledProcessError:
        print("Error metagene script")
        print("Moving on --")
        print("\033[91m{0}FAILURE ABOVE{0}\033[0m".format("-" * 30))


def find_clip_bed_files(folder):
    """
    Get the list of clip bed in a folder.

    :param folder: (str) a folder containing clip bed
    :return: (list of str) the list of clip bed in the folder.
    """
    cmd = "find %s -name *.bed.gz -type f" % folder
    result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    result = result.decode("ascii").split("\n")[:-1]
    return result


def find_wanted_exon_list(folder, sfname):
    """
    Get the up and down regulated exons files by ``sfname`` in ``folder``.

    :param folder: (str) a directory containing exon file
    :param sfname: (str) the name of a splicing factor
    :return: (list of 2 values) the list of wanted exon list.
    """
    cmd = "find %s -name input*%s*%s.txt"
    cmd_up = cmd % (folder, sfname, "up")
    cmd_down = cmd % (folder, sfname, "down")
    res_up = subprocess.check_output(cmd_up, shell=True,
                                     stderr=subprocess.STDOUT).decode("ascii")
    if len(res_up.replace("\n", "")) == 0:
        return (None, None), (None, None)
    res_up = res_up.split("\n")[:-1]
    res_down = subprocess.check_output(cmd_down, shell=True,
                                       stderr=subprocess.STDOUT)
    res_down = res_down.decode("ascii").split("\n")[:-1]
    res = res_up + res_down
    if len(res) != 2:
        raise IndexError("Error the len of res if %s while it should be 2" %
                         len(res))
    names = ("%s_up" % sfname, "%s_down" % sfname)
    return res, names


def find_wanted_and_ctrl_exon_list(folder, sfname):
    """
    Get the down regulated exons and every down-regulated exons \
     by at least one splicing factor files by ``sfname`` in ``folder``.

    :param folder: (str) a directory containing exon file
    :param sfname: (str) the name of a splicing factor
    :return: (list of 2 values) the list of wanted exon list.
    """
    cmd = "find %s -name input*%s*%s.txt"
    cmd_ctrl = "find %s -name CCE_exons.txt" % folder
    cmd_down = cmd % (folder, sfname, "down")
    res_ctrl = subprocess.check_output(
        cmd_ctrl, shell=True, stderr=subprocess.STDOUT).decode("ascii")
    if len(res_ctrl.replace("\n", "")) == 0:
        return (None, None), (None, None)
    res_ctrl = res_ctrl.split("\n")[:-1]
    res_down = subprocess.check_output(cmd_down, shell=True,
                                       stderr=subprocess.STDOUT)
    res_down = res_down.decode("ascii").split("\n")[:-1]
    res = res_ctrl + res_down
    if len(res) != 2:
        raise IndexError("Error the len of res if %s while it should be 2" %
                         len(res))
    names = ("CCE_exons", "%s_down" % sfname)
    return res, names


@lp.parse(input_folder="dir", folder_exon="dir", exon_bed="file",
          chrom_size_file="file", output="dir", metagene_script="file",
          ctrl=["y", "n"])
def main(input_folder, folder_exon, exon_bed, chrom_size_file, output,
         metagene_script, ctrl="n"):
    """

    :param input_folder: (str) a folder containing clip bed files files
    :param folder_exon: (str) a folder containing exon list
    :param exon_bed: (str) a bed file corresponding to fasterdb exons
    :param chrom_size_file: (str) a file indicating the length of every \
    human chromosome
    :param output: (str) path where the result will be created
    :param metagene_script: (str) path to the metagene script.
    :param ctrl: (str) 'y' to see only down exons and control exons \
    'n' to see up and down exons list
    """
    clip_beds = sorted(find_clip_bed_files(input_folder))
    for clip_bed in clip_beds:
        sfname = os.path.basename(clip_bed).split("_")[0].upper()
        print("Working on %s : files => %s" % (sfname, clip_bed))
        if ctrl == "n":
            exon_files, exonf_names = find_wanted_exon_list(folder_exon,
                                                            sfname)
        else:
            exon_files, exonf_names = find_wanted_and_ctrl_exon_list(
                folder_exon, sfname)
        if exon_files[0] is not None:
            color = "FF0000,1a8cff"
            for pos in ["n", "n-1", "n+1"]:
                coverage_wrapper(clip_bed, exon_bed, exon_files, exonf_names,
                                 chrom_size_file, output, metagene_script, "",
                                 color, "", pos)
        else:
            print("\tWarning : %s not found in exons folder" % sfname)


if __name__ == "__main__":
    main()
