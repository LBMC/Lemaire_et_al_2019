#!/usr/bin/python3

#####################
## Script to retrieve sequences for a list of annotations
## Usage: [python3] list_coord_seq_retriever.py <fasta_dir> <annot_file_list> <out_file_list>
# annot_file_list: list of files with annotations for which to retrieve the sequences ( ex: annot1,annot2,... )
#   ( annot_file format: "id ref:start-end strand" tab-separated, start is first base of the annotation )
# out_file_list: list of the output files in which to record the sequences of the corresponding annot_file ( ex: out1,out2,... )
#   ( annot1 > out1, annot2 > out2, ... )
#####################

import os,sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## parse the arguments
fasta_dir=sys.argv[1] #~ '/home/sebastien/analyses_SEBASTIEN/data/genome/hg19/EnsEMBL_all'
annot_file_list=sys.argv[2].split(',')
out_file_list=sys.argv[3].split(',')

## load the sequeces from the fasta files, one per chromosome (avoid the additional contigs ), from the given directory
# list the fasta files
for xxx in os.walk( fasta_dir):
    if xxx[ 0 ] == fasta_dir:
        files = [ yyy for yyy in xxx[ 2 ] if ( yyy[ -3: ] == '.fa' and yyy[ 0:2 ] != 'GL' ) ]
        break
        pass
    #~ print( xxx )
    pass

# load the sequence from each file
fasta_dict= {}
for xxx in files:
    fasta_dict[ xxx.split('.')[ 0 ] ] = SeqIO.read(fasta_dir + "/" + xxx, "fasta")
    pass


## retrieve the sequence for all the annotations and print to stdout
for index, annot_file_path in enumerate(annot_file_list):
    if annot_file_path == '-' and len( annot_file_list ) == 1:
        annot_file = sys.stdin.readlines()[ 1: ]
    else:
        annot_file = open( annot_file_path, 'r' ).readlines()[ 1: ]
        pass

    out_file_path = out_file_list[ index ]
    if out_file_path == '-' and len( out_file_list ) == 1:
        out_file = sys.stderr
    elif out_file_path == 'stdout':
        out_file = sys.stdout
    else:
        out_file = open( out_file_path, 'w' )
        pass

    print( 'coordinates\tstrand\tsequence\texon_id', file=out_file )
    for line in annot_file:
        exon_id, coord, strand = line.split()
        ref, borders = coord.split(':')
        start, end = [ int( xxx ) for xxx in borders.split('-') ]

        annot_seq = fasta_dict[ ref ].seq[ start-1:end ]
        if strand == '-1':
            annot_seq = annot_seq.reverse_complement()
            pass

        print( coord + '\t' + strand + '\t' + str( annot_seq ) + '\t' + exon_id, file=out_file )
        # if exon_id == '104890':
        #     break
        pass

##########
