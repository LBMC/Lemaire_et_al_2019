#!/usr/bin/python3

## Usage: python3 <script> <bw1,bw2,...> <bed1,bed2,...> <out_file_list> [ --add-chr 0,1,... ] [ --remove-duplicates ] [ --with-ids ]

import sys, os
import pyBigWig
from math import isnan
from statistics import mean
from retrieve_bw_values_binned_functions import *


####

bw_path_list = sys.argv[ 1 ].split( ',' )
bed_path_list = sys.argv[ 2 ].split( ',' )
out_file_list = sys.argv[ 3 ].split( ',' )

flag = '--add-chr'
if flag in sys.argv:
    add_chr_idx = sys.argv.index( flag )
    add_chr = [ int( xxx ) for xxx in sys.argv[ add_chr_idx + 1 ].split( ',' ) ]
    pass

nbins = None
flag = '--various-length'
if flag in sys.argv:
    nbins_idx = sys.argv.index( flag )
    nbins = int( sys.argv[ nbins_idx + 1 ] )
    pass

group_concat = False
flag = '--gene-concat'
if flag in sys.argv:
    group_concat = True
    pass

side_cut = [ 0, 0 ]
flag = '--side-cut'
if flag in sys.argv:
    side_cut = sys.argv[ sys.argv.index( flag ) + 1 ].split( ',' )
    side_cut = [ int( xxx ) - 1 for xxx in side_cut ]
    pass

## load bigwigs
bw_list = [ pyBigWig.open( bw_path ) for bw_path in bw_path_list ]

## load beds
bed_list = [ load_table( bed_path ) for bed_path in bed_path_list ]

## filter out duplicates
flag = '--remove-duplicates'
if flag in sys.argv:
    print( 'Remove duplicates !', file=sys.stderr )
    for bed_idx, bed in enumerate( bed_list ):
        print( 'BED idx: ' + str( bed_idx ), file=sys.stderr )
        print( '# line in bed: ' + str( len( bed ) ), file=sys.stderr )
        check_set = set()
        new_bed = []
        for line in bed:
            line_coord = ( line[ 0 ], line[ 1 ], line[ 2 ], line[ 5 ] )
            if line_coord not in check_set: 
                new_bed.append( line )
                check_set.add( line_coord )
                pass
            pass
        bed_list[ bed_idx ] = new_bed
        print( '# line in new bed: ' + str( len( bed_list[ bed_idx ] ) ), file=sys.stderr )
        pass
    pass


## get values
for bw_idx, bw in enumerate( bw_list ):
    for bed_idx, bed in enumerate( bed_list ):
        ## check for annotation of null length
        for line_idx, line in enumerate( bed ):
            if int( line[ 2 ] ) - int( line[ 1 ] ) == 0:
                print( 'Null length. ' + '\t'.join( line ) + '\nSkip it!', file=sys.stderr )
                #bed[ line_idx ] = None
                pass
        print( len( bed ), file=sys.stderr )
        bed = [ line for line in bed if int( line[ 2 ] ) - int( line[ 1 ] ) != 0 ]
        print( len( bed ), file=sys.stderr )

        ## get output file path
        if '--stdout' in sys.argv:
            out_file = sys.stdout
        else:
            out_file = open( out_file_list[ bw_idx * len( bed_list ) + bed_idx ], 'w' )
            pass

        ## check in bw with ref in 'chr' format
        if '--add-chr' in sys.argv:
            bed = check_ref_chr( bed, add_chr[ 0 ] )
            pass

        ## retrieve values
        #~ all_means = retrieve_multi_mean_cov( bed, bw, nbins=nbins, group_concat=group_concat )
        all_means = retrieve_multi_interval( bed, bw, nbins=nbins, group_concat=group_concat, start_left=side_cut[ 0 ], end_right=side_cut[ 1 ], mean_val=True )
        flag = '--with-ids'
        if flag in sys.argv:
            if len( bed[ 0 ] ) >= 4:
                annot_id_list = [ xxx[ 3 ] for xxx in bed ]
            else:
                annot_id_list = [ str( xxx + 1 ) for xxx in range( len( bed ) ) ]
                pass

            for mean_idx, mean_val in enumerate( all_means ):
                print( '\t'.join( [ annot_id_list[ mean_idx ], str( mean_val ) ] ), file=out_file )
                pass
        else:
            print( '\t'.join( [ str( xxx ) for xxx in all_means ] ), file=out_file )
            pass

        if not '---stdout' in sys.argv:
            out_file.close()
            pass
        pass
    pass



####
