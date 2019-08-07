#!/usr/bin/env python3


"""
Description:
    The goal of this script is to create a figure showing the coverage of a bed file in the form of \
    a metaexon graphics.

    Requirements: betools and bedGraphToBigWig must be in your path.
"""

import subprocess
import os
import argparse
import union_dataset_function
import sqlite3
import gzip
import group_factor
import re


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
    bed_chrom = outfolder + os.path.basename(bed_file).replace(".bed", "nochr.bed").replace(".gz", "")
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
    cmd_bg = "bedtools genomecov -bg -i %s -g %s > %s" % (sort_bed, refsize, bg_file)
    print(cmd_bg)
    subprocess.check_call(cmd_bg, shell=True, stderr=subprocess.STDOUT)
    # transformation to bigwig
    cmd_bw = "bedGraphToBigWig %s %s %s" % (bg_file, refsize, bw_file)
    print(cmd_bw)
    subprocess.check_call(cmd_bw, shell=True, stderr=subprocess.STDOUT)
    return bw_file


def metagene_covergae_launcher(metagene_script, bw_file, cond, bed_annot, annot_name, color, outname, outdir, off_set):
    """
    Launch the metagene script program.

    :param metagene_script: (string) path where the metagene script is located
    :param bw_file: (string) path to the bigwigfile
    :param cond: (string) a title
    :param bed_annot: (string) a bed file containing the chosen annotation to work on for metagene creation
    :param annot_name: (string) the name of the annotation file
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    :param outdir:  (string) path where the result will be created
    :param off_set: (int) offset for 0 in the graphic

    """
    cmd = """Rscript %s -bw %s -cond %s -rep 1 -annot %s -prefixes %s --color-pallette %s -sign %s\
    -out_dir %s -off_set %s -comp_pair %s""" % (metagene_script, bw_file, cond, bed_annot, annot_name, color,
                                                outname, outdir, off_set, cond)
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


def get_control_exon_information(cnx, exon_type):
    """
    Get the gene symbol, the gene id and the position of every ``exon_type`` exons in fasterDB.

    :param cnx: (sqlite3 object) allow connection to sed database
    :param exon_type: (string) the type of control exon we want to use
    :return: (list of 2 int) gene_id and exon_pos
    """
    cursor = cnx.cursor()
    if exon_type != "ALL":
        query = """SELECT gene_id, exon_pos
                   FROM sed
                   WHERE exon_type LIKE '%{}%'
                   """.format(exon_type)
    else:
        query = """SELECT  gene_id, exon_pos
                   FROM sed
                   AND t1.id_gene = t2.id
                """
    cursor.execute(query)
    result = cursor.fetchall()
    # turn tuple into list
    nresult = []
    for exon in result:
        nresult.append(list(exon))
    return nresult


def get_gene_id(cnx, coord, gene_symbol):
    """
    From coordinates of an exon retrieve the gene_id of the gene that contains this exons.

    :param cnx: (pymysql connect instance) connection to splicing lore database
    :param coord: (string) the coordinates of an exon (format chr:start-stop)
    :param gene_symbol: (string) the gene name of the exon with coordinates ``coord``
    :return: (int) the gene id containing the exon with the coordinates ``coord``
    """
    cursor = cnx.cursor()
    coords = re.split(":|-", coord)
    query = """SELECT t1.id_gene
               FROM exons_genomiques_bis t1, genes t2
               WHERE t1.id_gene = t2.id
               AND t2.official_symbol = '%s'
               AND t1.chromosome = '%s'
               AND t1.start_sur_chromosome = %s
               AND t1.end_sur_chromosome = %s""" % (gene_symbol, coords[0], coords[1], coords[2])
    cursor.execute(query)
    result = cursor.fetchall()
    if len(result) > 1:
        print("warning : more than one exon were found with the coordinates %s" % coord)
    if not result[0]:
        print("error :  no exon was found with the coordinates %s" % coord)
        exit(1)
    return result[0][0]


def get_exon_list(cnx, annotation_name, regulation):
    """
    Get the exon_list wanted.

    :param cnx: (sqlite3 connect object) connection to sed database
    :param annotation_name: (string) GC-AT or a sf_name
    :param regulation: (string) the regulation of an exon list by a factor(s)
    :return: (list of 2 int) gene id and exon_pos
    """
    if "GC" in annotation_name or "AT" in annotation_name:
        annotation_name = annotation_name.split("_")[0]
        folder = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/")
        my_file = "%s%s_rich_exons" % (folder, annotation_name)
        exon_list = extract_exon_list(my_file)
    elif "U1-FACTORS" in annotation_name or "U2-FACTORS" in annotation_name:
        annotation_name = annotation_name.split("_")[0]
        dic_name = {"U1-FACTORS": ["SNRPC", "SNRNP70", "DDX5_DDX17"], "U2-FACTORS": ["U2AF2", "SF1", "SF3A3", "SF3B4"]}
        exon_list = []
        for sf_name in dic_name[annotation_name]:
            exon_list += union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
        exon_list = union_dataset_function.washing_events_all(exon_list)
    else:
        annotation_name = annotation_name.split("_")[0]
        sf_name = annotation_name.upper()
        sf_name = sf_name.replace("SFRS", "SRSF").replace("TRA2A", "TRA2A_B").replace("DDX5-17", "DDX5_DDX17")
        exon_list = union_dataset_function.get_every_events_4_a_sl(cnx, sf_name, regulation)
    return exon_list


def handle_bed_creator(cnx, cnx_fasterdb, annotation_name, template_fodler, chrom_size_file, regulation):
    """
    Get the name of the annotation file.

    :param cnx: (sqlite3 connect object) connection to sed database
    :param cnx_fasterdb: (sqtite3 connect object) connection to fasterDB
    :param annotation_name: (string) GC-AT or a sf_name
    :param template_fodler: (string) a folder containing template bed file.
    :param chrom_size_file: (string) a file containing the length of hg19 chromosomes.
    :param regulation: (string) up or down
    :return: 4 elements:

        1. (string) the list of bed files (coma separated) corresponding to intron-exon junctions containing the 3'SS
        2. (string) the list of bed files (coma separated) corresponding to intron-exon junctions containing the 5'SS
        3. (string) the list of name (coma separated) for the bed files ) corresponding to intron-exon junctions \
        containing the 3'SS
        4. (string) the list of name (coma separated) for the bed files ) corresponding to intron-exon junctions \
        containing the 5'SS
    """
    final_template_start = []
    final_template_stop = []
    if annotation_name == "GC-AT":
        names = ["GC_group", "AT_group"]
    else:
        sf_list = annotation_name.split(",")
        for i in range(len(sf_list)):
            sf_list[i] = sf_list[i].upper()
            sf_list[i] = sf_list[i].replace("SFRS", "SRSF").replace("TRA2A", "TRA2A_B")
        names = ["%s_%s" % (sf_name, regulation) for sf_name in sf_list]
    for i in range(len(names)):
        if not os.path.isfile("%s%s_%s_add_i50_o200.bed" % (template_fodler, names[i], "start")):
            exon_list = get_exon_list(cnx, names[i], regulation)
            templates_bed = bed_creator(cnx_fasterdb, exon_list, template_fodler, names[i], chrom_size_file)
        else:
            start = "%s%s_%s_add_i50_o200.bed" % (template_fodler, names[i], "start")
            templates_bed = [start, start.replace("start", "stop")]
            print("recovering template files : %s" % start)
        final_template_start.append(templates_bed[0])
        final_template_stop.append(templates_bed[1])
    finale_names_start = []
    finale_names_stop = []
    for my_name in names:
        finale_names_start.append(my_name + "_start")
        finale_names_stop.append(my_name + "_stop")
    finale_names_start = ",".join(finale_names_start)
    finale_names_stop = ",".join(finale_names_stop)
    final_templates_start = ",".join(final_template_start)
    final_templates_stop = ",".join(final_template_stop)
    return final_templates_start, final_templates_stop, finale_names_start, finale_names_stop


def get_exon_list_file(list_file):
    """

    :param list_file: (str) a file containing exons
    :return: (list of 2 int) a list of exons
    """
    exon_list = []
    with open(list_file, "r") as infile:
        for line in infile:
            line = line.replace("\n", "")
            line = line.split("\t")
            exon_list.append(list(map(int, line)))
    return exon_list


def handle_bed_creator_file(cnx_fasterdb, list_file, name_file,
                            template_fodler, chrom_size_file):
    """
    Get the name of the annotation file.

    :param cnx_fasterdb: (sqtite3 connect object) connection to fasterDB
    :param list_file: (list of str) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (list of str) the name of each files of exons \
    given in ``list_file``
    :param template_fodler: (string) a folder containing template bed file.
    :param chrom_size_file: (string) a file containing the length of hg19 chromosomes.
    :param regulation: (string) up or down
    :return: 4 elements:

        1. (string) the list of bed files (coma separated) corresponding to intron-exon junctions containing the 3'SS
        2. (string) the list of bed files (coma separated) corresponding to intron-exon junctions containing the 5'SS
        3. (string) the list of name (coma separated) for the bed files ) corresponding to intron-exon junctions \
        containing the 3'SS
        4. (string) the list of name (coma separated) for the bed files ) corresponding to intron-exon junctions \
        containing the 5'SS
    """
    final_template_start = []
    final_template_stop = []
    for i in range(len(name_file)):
        if not os.path.isfile("%s%s_%s_add_i50_o200.bed" % (template_fodler, name_file[i], "start")):
            exon_list = get_exon_list_file(list_file[i])
            templates_bed = bed_creator(cnx_fasterdb, exon_list, template_fodler, name_file[i], chrom_size_file)
        else:
            start = "%s%s_%s_add_i50_o200.bed" % (template_fodler, name_file[i], "start")
            templates_bed = [start, start.replace("start", "stop")]
            print("recovering template files : %s" % start)
        final_template_start.append(templates_bed[0])
        final_template_stop.append(templates_bed[1])
    finale_names_start = []
    finale_names_stop = []
    for my_name in name_file:
        finale_names_start.append(my_name + "_start")
        finale_names_stop.append(my_name + "_stop")
    finale_names_start = ",".join(finale_names_start)
    finale_names_stop = ",".join(finale_names_stop)
    final_templates_start = ",".join(final_template_start)
    final_templates_stop = ",".join(final_template_stop)
    return final_templates_start, final_templates_stop, finale_names_start, finale_names_stop

def load_chrom_size_dic(chrom_size_file):
    """
    Load in the form a a dictionary the chromosome size file.

    :param chrom_size_file: (float) chromosome size filz
    :return: (dictionary of int) dictionary that links each chromosome to it's size
    """
    dic_size = {}
    with open(chrom_size_file, "r") as size_file:
        for line in size_file:
            line = line.replace("\n", "").split("\t")
            if "_" not in line[0]:
                dic_size[line[0]] = int(line[1])
    return dic_size


def bed_creator(cnx_fasterdb, exon_list, dest_folder, name_bed, chrom_size_file):
    """
    Create a bed file containing all the exons regulated by ``sf_name`` with ``regulation``.

    :param cnx_fasterdb: (sqlite3 connector object) connection to fasterDB lite database
    :param exon_list: (list of 2 int) gene id + exon pos
    :param dest_folder: (string) path where the bed will be created
    :param name_bed: (string) name of the file
    :param chrom_size_file: (string) a file containing chromosome size
    :return: (string) the name of the bed file created
    """
    in_exon = 50 - 1
    out_exon = 200
    dic_size = load_chrom_size_dic(chrom_size_file)
    cursor = cnx_fasterdb.cursor()
    exon_info_left = []
    exon_info_right = []
    for exon in exon_list:
        query = """SELECT t1.chromosome, t1.start_on_chromosome, t1.end_on_chromosome, 
                          t2.official_symbol, t1.pos_on_gene,
                   t2.strand, t1.end_on_chromosome - t1.start_on_chromosome + 1
                   FROM exons t1, genes t2
                   WHERE t1.id_gene = t2.id
                   AND t1.id_gene = %s
                   AND t1.pos_on_gene = %s
                   ORDER BY t1.chromosome ASC, t1.start_on_chromosome ASC
                """ % (exon[0], exon[1])
        cursor.execute(query)
        res = cursor.fetchall()
        if len(res) > 1:
            print("Error, only one exon should be found for %s_%s exon" % (exon[0], exon[1]))
            exit(1)
        exon_left = None
        exon_right = None
        if res[0][5] == 1:
            if res[0][1] - out_exon > 0 and res[0][1] + in_exon <= dic_size[str(res[0][0])]:
                exon_left = [res[0][0], res[0][1] - out_exon, res[0][1] + in_exon] +\
                            ["%s_%s" % (res[0][3], res[0][4])] + ["."] + [res[0][5], res[0][6]]
            if res[0][2] - in_exon > 0 and res[0][2] + out_exon <= dic_size[str(res[0][0])]:
                exon_right = [res[0][0], res[0][2] - in_exon, res[0][2] + out_exon] + \
                             ["%s_%s" % (res[0][3], res[0][4])] + ["."] + [res[0][5], res[0][6]]
        else:
            if res[0][2] - in_exon > 0 and res[0][2] + out_exon <= dic_size[str(res[0][0])]:
                exon_left = [res[0][0], res[0][2] - in_exon, res[0][2] + out_exon] + \
                            ["%s_%s" % (res[0][3], res[0][4])] + ["."] + [res[0][5], res[0][6]]
            if res[0][1] - out_exon > 0 and res[0][1] + in_exon <= dic_size[str(res[0][0])]:
                exon_right = [res[0][0], res[0][1] - out_exon, res[0][1] + in_exon] + \
                             ["%s_%s" % (res[0][3], res[0][4])] + ["."] + [res[0][5], res[0][6]]
        if exon_left:
            exon_info_left.append("\t".join(list(map(str, exon_left))))
        if exon_right:
            exon_info_right.append("\t".join(list(map(str, exon_right))))
    contents = ["\n".join(exon_info_left), "\n".join(exon_info_right)]
    names = ["start", "stop"]
    final_names = []
    for i in range(len(contents)):
        filename = "%s%s_%s_add_i50_o200.bed" % (dest_folder, name_bed, names[i])
        # final_name = filename.replace(".bed", "_add_i50_o200.bed")
        final_names.append(filename)
        with open(filename, "w") as bedfile:
            bedfile.write(contents[i] + "\n")
        # if names[i] == "start":
        #     fasterdb_bed_add_exon_type.add_intron_sequence(filename, final_name, chrom_size_file, 200, 50)
        # else:
        #     fasterdb_bed_add_exon_type.add_intron_sequence(filename, final_name, chrom_size_file, 50, 200)
    return final_names


def main(input_bed, refsize, output, metagene_script, cond, annotation, color, outname, regulation):
    """

    :param input_bed: (string) a bed file
    :param refsize: (string) a file containing the hg19 chromosome size.
    :param output: (string) path where the result will be created.
    :param metagene_script: (string) path where the metagene script is located
    :param annotation: (string) the annotation wanted
    :param cond: (string) a title
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    :param regulation: (string) up or down
    """
    seddb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/sed.db")
    fasterdb = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/fasterDB_lite.db")
    template_fodler = seddb.replace("data/sed.db", "result/template/")
    if not os.path.isdir(template_fodler):
        os.mkdir(template_fodler)
    cnx = sqlite3.connect(seddb)
    cnx_fasterdb = sqlite3.connect(fasterdb)
    tplt_start, tplt_stop, n_start, n_stop = handle_bed_creator(cnx, cnx_fasterdb, annotation,
                                                                template_fodler, refsize, regulation)
    output_bw = output + "bigwig/"
    output_graph = output + "figure/"
    if not os.path.isdir(output_bw):
        os.mkdir(output_bw)
    if not os.path.isdir(output_graph):
        os.mkdir(output_graph)
    if not color:
        prefixe = n_start.split(",")
        color = [group_factor.color_dic[my_name.split("_")[0].replace("retained", "CCE")] for my_name in prefixe]
        color = ",".join(color).replace("#", "")

    bw_file = bed2bw(input_bed, refsize, output_bw)
    if not outname:
        outname = ""
        if annotation == "GC-AT":
            outname += "GC_AT_group_down"
        else:
            outname += "spliceosome_%s" % regulation
    try:
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_start,
                                   n_start, color, outname + "_start", output_graph, 200)
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_stop, n_stop,
                                   color, outname + "_stop", output_graph, 50)
    except subprocess.CalledProcessError:
        print("Error metagene script")


def figure_2h(input_bed, list_file, name_file, refsize, seddb, fasterdb, output,
              metagene_script, cond, color, outname):
    """
    Create the figure 2f.

    :param input_bed: (string) a bed file
    :param list_file: (list of str) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (list of str) the name of each files of exons \
    given in ``list_file``
    :param refsize: (string) a file containing the hg19 chromosome size.
    :param output: (string) path where the result will be created.
    :param metagene_script: (string) path where the metagene script is located
    :param cond: (string) a title
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    """
    template_fodler = seddb.replace("data/sed.db", "result/template/")
    if not os.path.isdir(template_fodler):
        os.mkdir(template_fodler)

    cnx = sqlite3.connect(seddb)
    cnx_fasterdb = sqlite3.connect(fasterdb)
    tplt_start, tplt_stop, n_start, n_stop = \
        handle_bed_creator_file(cnx_fasterdb, list_file, name_file,
                                template_fodler, refsize)
    output_bw = output + "bigwig/"
    output_graph = output + "figure/"
    if not os.path.isdir(output_bw):
        os.mkdir(output_bw)
    if not os.path.isdir(output_graph):
        os.mkdir(output_graph)
    if not color:
        color_dic = group_factor.color_dic
        prefixe = n_start.split(",")
        try:
            color = [color_dic[my_name.split("_")[0].replace("retained", "CCE")] for my_name in prefixe]
        except KeyError:
            color = [color_dic["GC_pure"], color_dic["AT_pure"], color_dic["CCE"]]
            color = color[:len(prefixe)]
        color = ",".join(color).replace("#", "")

    bw_file = bed2bw(input_bed, refsize, output_bw)
    if not outname:
        outname = "_".join(name_file)
    try:
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_start,
                                   n_start, color, outname + "_start", output_graph, 200)
        metagene_covergae_launcher(metagene_script, bw_file, cond, tplt_stop, n_stop,
                                   color, outname + "_stop", output_graph, 50)
    except subprocess.CalledProcessError:
        print("Error metagene script")
    cnx.close()
    cnx_fasterdb.close()


def main_2h(bed_folder, list_file, name_file, refsize, seddb, fasterdb, output,
            metagene_script, cond, color, outname, num_fig):
    """

    :param bed_folder: (string) a folder containing bed files
    :param list_file: (list of str) list of exons files in the form \
    of GC_rich_exon file.
    :param name_file: (list of str) the name of each files of exons \
    given in ``list_file``
    :param refsize: (string) a file containing the hg19 chromosome size.
    :param output: (string) path where the result will be created.
    :param metagene_script: (string) path where the metagene script is located
    :param cond: (string) a title
    :param color: (string) colors chosen
    :param outname: (string) path of the output file
    :param num_fig: (str) the num of the fig
    """

    bed_list = os.listdir(bed_folder)
    for input_bed in bed_list:
        print("Working on %s" % input_bed)
        noutput = output + "/" + num_fig + "_" + \
            input_bed.replace(".bed.gz", "") + "/"
        if not os.path.isdir(noutput):
            os.mkdir(noutput)
        figure_2h("%s/%s" % (bed_folder, input_bed), list_file, name_file,
                  refsize, seddb, fasterdb, noutput, metagene_script, cond,
                  color, outname)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    This program will convert a bed file to a BigWig file and will launch the metagene coverage program
    """)
    # Arguments for the parser

    req_arg = parser.add_argument_group("Required arguments")

    req_arg.add_argument('--input', dest='input', help="An input bed file for which we want to produce an output graph",
                         required=True)
    req_arg.add_argument('--refsize', dest='refsize',  help="The size of chromosome of hg19", required=True)
    parser.add_argument('--output', dest='output', help="path where the result will be created", default="./")
    req_arg.add_argument('--metagene_script', dest='metagene_script', help="path where the meatagne script is located",
                         required=True)
    parser.add_argument('--title', dest='title', help="title of the graphic",
                        default=None)
    req_arg.add_argument('--annotation', dest='annotation',
                         help="GC-AT or the name of a splicing factor", required=True)
    parser.add_argument('--color', dest='color', help="color of the annotation", default=None)
    parser.add_argument('--outname', dest='outname', help="the name of the figure", default=None)
    parser.add_argument('--regulation', dest='regulation', help="regulation", default="down")
    args = parser.parse_args()

    if os.path.isdir(args.input):
        list_files = os.listdir(args.input)
        for my_file in list_files:
            print("Working on %s" % my_file)
            if not os.path.isdir(args.output):
                os.mkdir(args.output)
            output = args.output + my_file.replace(".bed.gz", "") + "/"
            if not os.path.isdir(output):
                os.mkdir(output)

            # main(args.input + my_file, args.refsize, output, args.metagene_script, args.title, "GC-AT",
            #  args.color, args.outname, args.regulation)
            # annot = "SNRPC,SNRNP70,U2AF2,SF1"
            main(args.input + my_file, args.refsize, output, args.metagene_script, args.title, args.annotation,
                 args.color, args.outname, "down")
            # main(args.input + my_file, args.refsize, output, args.metagene_script, args.title,
            #      args.annotation,
            #      args.color, args.outname, "up")
    else:
        main(args.input, args.refsize, args.output, args.metagene_script, args.title, args.annotation,
             args.color, args.outname, args.regulation)


if __name__ == "__main__":
    launcher()
