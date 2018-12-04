#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import get_bed_sequences
import bed_handler
import argparse
import os


def main(clip_bed, bed_folder, hg19_reference, output, size, prg_path, enrichment):
    """
    Create a weblogo of the sequences in ``clip_bed`` overlapping those in ``fasterdb_bed``.

    :param clip_bed: (string) a bed containing motifs bing by a particular splicing factor
    :param bed_folder: (string)a folder containing interest bed files
    :param hg19_reference: (string) a fasta file containing the sequences of every chromosome in hg19 or \
    a dictionary containing the sequences of hg19
    :param output: (string) path where the weblogo will be created
    :param size: (int) the size of the weblogo
    :param prg_path: (string) path to motif finder programs
    :param enrichment: (string) y/n : y to perform en enrichment motif analysis and N to \
    produce only a weblogo centered in the clip sequence
    """
    if isinstance(hg19_reference, str):
        print("Loading hg19 genome !")
        dic_records = get_bed_sequences.create_sequence_dic(hg19_reference)
    else:
        print(" --- ANALYSING %s ---" % clip_bed)
        dic_records = hg19_reference
    outputs = [output + "1based_bed/", output + "intesect_bed/", output + "frequencies/"]
    for cur_output in outputs:
        if not os.path.isdir(cur_output):
            os.mkdir(cur_output)
    print("Creation of a one based bed")
    bed_list = [bed_handler.bed_0based_2_1based(clip_bed, outputs[0])]
    print("Creation of intersect bed")
    list_template = bed_handler.get_files(bed_folder)
    for bed in list_template:
        base_name = os.path.basename(bed).split(".")[0]
        bed_list.append(bed_handler.intersect_bed(bed, bed_list[0], outputs[1], base_name))
    for i in range(len(bed_list)):
        print("Working on %s bed file" % bed_list[i])
        list_coordinates = get_bed_sequences.get_bed_coordinates(bed_list[i], dic_records.keys())
        list_coordinates = list_coordinates
        print("   --> nb list_coordinates: %s" % len(list_coordinates))
        list_sequence = get_bed_sequences.get_middle_sequences(list_coordinates, dic_records, size)
        dic_freq = get_bed_sequences.get_nt_frequencies(list_sequence)
        name_result = os.path.basename(bed_list[i]).split(".")[1]
        with open("%sfrequencies_%s.py" % (outputs[2], name_result), "w") as outfile:
            outfile.write("dic_freq = %s\n" % str(dic_freq))
        if enrichment == "Y":
            motif_search_dir = output + "motif_search_%s/" % name_result
            if not os.path.isdir(motif_search_dir):
                os.mkdir(motif_search_dir)
            fa_file = get_bed_sequences.write_fasta_file(list_coordinates, dic_records, motif_search_dir)
            print("Searching enriched motifs")
            bed_handler.meme_launcher(prg_path, fa_file, motif_search_dir)
            zagros_size = min(size, 8)
            bed_handler.zagros_launcher(prg_path, fa_file, zagros_size, motif_search_dir)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    This program aims to create a weblogo from the feature in a bed file obtained with a clip experiment \
    if those features overlap another bedfile corresponding to the fasterdb exons with 200 additional bp at \
    each side of the exon.
    """)
    # Arguments for the parser

    req_arg = parser.add_argument_group("Required arguments")

    req_arg.add_argument('--clip_bed', dest='clip_bed', help="bed coming from a clip experiment or a folder containing "
                                                             "a lot of bed",
                         required=True)
    parser.add_argument('--enrichment', dest='enrichment',
                        help="(Y/N)Y perform a motif search in clipped sequence (time consuming), "
                             "N if you want only a belogo centered in the clip sequence", default="N")
    parser.add_argument('--fasterdb_bed', dest='fasterdb_bed', help="folder where are located the bed'",
                        default=os.path.realpath(os.path.dirname(__file__)).
                        replace("src", "data/bed_template/"))
    parser.add_argument('--hg19_reference', dest='hg19_reference',
                        help="fasta file containing the sequences of every chromosome in hg19",
                        default=os.path.realpath(
                            os.path.dirname(__file__)).replace("src",
                                                               "data/Homo_sapiens.GRCh37.dna.primary_assembly.fa"))
    req_arg.add_argument('--output', dest='output', help="path where the weblogo will be created",
                         required=True)
    parser.add_argument('--size', dest='size', help='size', default=10)
    parser.add_argument('--nb_iteration', dest='nb_iteration', help='nb control iteration', default=1000)
    parser.add_argument('--exon_type', dest='exon_type', help='exon type (ACE, CCE, or ALL)', default="CCE")
    parser.add_argument('--meme_path', dest="meme_path", help="path to meme program",
                        default=os.path.realpath(os.path.dirname(__file__)).replace("src", "data/meme_program/"))
    parser.add_argument('--chr_size', dest="chr_size", help="a tab file indicating the size of each chromosome in hg19",
                        default=os.path.realpath(os.path.dirname(__file__)).replace("src", "data/hg19.ren.chrom.sizes"))

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        parser.error("The output directory doesn't exist !")
    if args.output[-1] != "/":
        args.output += "/"

    if args.meme_path[-1] != "/":
        args.meme_path += "/"

    if not os.path.isfile(args.clip_bed) and not os.path.isdir(args.clip_bed):
        parser.error("The clip bed doesn't exist !")

    try:
        args.size = int(args.size)
        if args.size < 0:
            parser.error("ERROR : wrong weblogo size")
    except ValueError:
        parser.error("ERROR : wrong weblogo size")

    response = {"yes": "Y", "y": "Y", "Y": "Y", "no": "N", "n": "N", "N": "N"}
    if args.enrichment not in response.keys():
        parser.error("Wrong value for --enrichment, the values Y/N/y/n/yes/no are accepted")
    else:
        args.enrichment = response[args.enrichment]

    try:
        args.nb_iteration = int(args.nb_iteration)
        if args.nb_iteration < 1:
            parser.error("Wrong parameter for nb_iteration, it must be a positive value")
    except ValueError:
        parser.error("Wrong parameter for nb_iteration, it must be a number")

    if os.path.isfile(args.clip_bed):
        main(args.clip_bed, args.fasterdb_bed, args.hg19_reference, args.output, args.size, args.meme_path,
             args.enrichment)
    else:
        print("Loading hg19 genome !")
        dic_records = get_bed_sequences.create_sequence_dic(args.hg19_reference)
        list_bed = bed_handler.get_files(args.clip_bed)
        for mybed in list_bed:
            output = args.output + os.path.basename(mybed).split(".")[0] + "/"
            if not os.path.isdir(output):
                os.mkdir(output)
            main(mybed, args.fasterdb_bed, dic_records, output, args.size, args.meme_path,
                 args.enrichment)


if __name__ == "__main__":
    launcher()
