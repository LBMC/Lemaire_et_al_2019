#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import get_bed_sequences
import figure_maker
import bed_handler
import argparse
import os
import create_control_sets2 as create_control_sets

#import create_control_sets


def main(clip_bed, fasterdb_bed, hg19_reference, output, size, prg_path, enrichment, nb_iteration, exon_type,
         hg19_chromosome_file):
    """
    Create a weblogo of the sequences in ``clip_bed`` overlapping those in ``fasterdb_bed``
    :param clip_bed: (string) a bed containing motifs bing by a particular splicing factor
    :param fasterdb_bed: (string) a bed containing every fasterDB exons with 200bp of interval
    :param hg19_reference: (string) a fasta file containing the sequences of every chromosome in hg19
    :param output: (string) path where the weblogo will be created
    :param size: (int) the size of the weblogo
    :param prg_path: (string) path to motif finder programs
    :param enrichment: (string) y/n : y to perform en enrichment motif analysis and N to \
    produce only a weblogo centered in the clip sequence
    :param nb_iteration: (int) the chosen number of iteration to generated the control sequence
    :param exon_type: (string) the control exon chosen (CCE, ACE or ALL)
    :param hg19_chromosome_file: (string) a file containing hg19 chromosome size
    """
    project_name = os.path.basename(clip_bed).split(".")[0]
    outputs =[output + "1based_bed/", output + "intesect_bed/", output + "figures/"]
    for cur_output in outputs:
        if not os.path.isdir(cur_output):
            os.mkdir(cur_output)
    print("Creation of a one based bed")
    b1_bed = bed_handler.bed_0based_2_1based(clip_bed, outputs[0])
    # print("intersection with fasterDB bed")
    # intersect_bed = bed_handler.intersect_bed(fasterdb_bed, b1_bed, outputs[1])
    print("Loading hg19 genome !")
    dic_records = get_bed_sequences.create_sequence_dic(hg19_reference)
    # print("getting weblogo sequences")
    # list_coordinates = get_bed_sequences.get_bed_coordinates(intersect_bed)
    list_coordinates = get_bed_sequences.get_bed_coordinates(b1_bed, dic_records.keys())
    print("nb list_coordinates")
    print(len(list_coordinates))
    # list_sequence = get_bed_sequences.get_middle_sequences(list_coordinates, dic_records, size)
    # print("Making weblogo sequence")
    # figure_maker.web_logo_creator(list_sequence, project_name, outputs[2])
    print("Creating control exons sets")
    create_control_sets.set_debug(0)
    mean_ctrl_dic = create_control_sets.make_control_sets(nb_iteration, hg19_chromosome_file, fasterdb_bed, exon_type,
                                          list_coordinates, dic_records)
    print("creating test_set")
    # print("coord test : %s" % str(list_coordinates[0:3]))
    # mean_test_dic = create_control_sets.create_mean_frequency_dic(create_control_sets.get_sequences_frequencies(list_coordinates[0:3], dic_records))
    mean_test_dic = create_control_sets.create_mean_frequency_dic(
        create_control_sets.get_sequences_frequencies_test(list_coordinates, dic_records))
    figure_maker.make_control_figures(mean_test_dic, mean_ctrl_dic, outputs[2], exon_type, project_name, nb_iteration)

    if enrichment == "Y":
        motif_search_dir = output + "motif_search/"
        if not os.path.isdir(motif_search_dir):
            os.mkdir(motif_search_dir)
        fa_file = get_bed_sequences.write_fasta_file(list_coordinates, dic_records, motif_search_dir)
        del (dic_records)
        print("Searching enriched motifs")
        bed_handler.meme_launcher(prg_path, fa_file, motif_search_dir)
        bed_handler.zagros_launcher(prg_path, fa_file, size, motif_search_dir)


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

    req_arg.add_argument('--clip_bed', dest='clip_bed', help="bed coming from a clip experiment",
                         required=True)
    parser.add_argument('--enrichment', dest='enrichment', help="(Y/N)Y perform a motif search in clipped sequence (time consuming), "
                                                                 "N if you want only a belogo centered in the clip sequence",
                        default="N")
    parser.add_argument('--fasterdb_bed', dest='fasterdb_bed', help="bed file containing the fasterdb exons with "
                                                                    "200 additional bp at each side of the exon.",
                        default=os.path.realpath(os.path.dirname(__file__)).replace("src",
                                                                                     "data/fasterDB_exons_add200nt.bed"))
    parser.add_argument('--hg19_reference', dest='hg19_reference',
                        help="fasta file containing the sequences of every chromosome in hg19",
                        default=os.path.realpath(
                            os.path.dirname(__file__)).replace("src","data/Homo_sapiens.GRCh37.dna.primary_assembly.fa"))
    req_arg.add_argument('--output', dest='output', help="path where the weblogo will be created",
                        required=True)
    parser.add_argument('--size', dest='size', help='weblogo size', default=6)
    parser.add_argument('--nb_iteration', dest='nb_iteration', help='nb control iteration', default=1000)
    parser.add_argument('--exon_type', dest='exon_type', help='exon type (ACE, CCE, or ALL)', default="CCE")
    parser.add_argument('--meme_path', dest="meme_path", help="path to meme program",
                        default=os.path.realpath(os.path.dirname(__file__)).replace("src","data/meme_program/"))
    parser.add_argument('--chr_size', dest="chr_size", help="a tab file indicating the size of each chromosome in hg19",
                        default=os.path.realpath(os.path.dirname(__file__)).replace("src", "data/hg19.ren.chrom.sizes"))

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        parser.error("The output directory doesn't exist !")
    if args.output[-1] != "/":
        args.output += "/"

    if args.meme_path[-1] != "/":
        args.meme_path += "/"

    if not os.path.isfile(args.clip_bed):
        parser.error("The clip bed doesn't exist !")

    try:
        args.size = int(args.size)
        if args.size < 0:
            parser.error("ERROR : wrong weblogo size")
    except:
        parser.error("ERROR : wrong weblogo size")

    response = {"yes": "Y", "y": "Y", "Y":"Y", "no": "N", "n": "N", "N": "N"}
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

    main(args.clip_bed, args.fasterdb_bed, args.hg19_reference, args.output, args.size, args.meme_path, args.enrichment,
         args.nb_iteration, args.exon_type, args.chr_size)


if __name__ == "__main__":
    launcher()

