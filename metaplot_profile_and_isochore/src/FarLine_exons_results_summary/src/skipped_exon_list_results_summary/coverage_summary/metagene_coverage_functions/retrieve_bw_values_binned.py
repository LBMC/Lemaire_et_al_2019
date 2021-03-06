#!/usr/bin/python3

## Usage: python3 <script> <bw1,bw2,...> <bed1,bed2,...> <out_dir> <nbins> [ --add-chr 0,1,... ]

import sys, os
import pyBigWig
from math import isnan
from statistics import mean
from retrieve_bw_values_binned_functions import *


####

bw_path_list = sys.argv[ 1 ].split( ',' )
bed_path_list = sys.argv[ 2 ].split( ',' )
out_file_list = sys.argv[ 3 ].split( ',' )
nbins = int( sys.argv[ 4 ] )

flag = '--add-chr'
if flag in sys.argv:
    add_chr_idx = sys.argv.index( flag )
    add_chr = [ int( xxx ) for xxx in sys.argv[ add_chr_idx + 1 ].split( ',' ) ]
    pass

group_concat = False
flag = '--group-concat'
if flag in sys.argv:
    group_concat = True
    pass

## load bigwigs
bw_list = [ pyBigWig.open( bw_path ) for bw_path in bw_path_list ]

## load beds
bed_list = [ load_table( bed_path ) for bed_path in bed_path_list ]


## get values
for bw_idx, bw in enumerate( bw_list ):
    for bed_idx, bed in enumerate( bed_list ):
        ## get output file path
        if '--stdout' in sys.argv:
            out_file = sys.stdout
        else:
            out_file = open( out_file_list[ bw_idx * len( bed_list ) + bed_idx ], 'w' )
            pass

        ## check for annotation of null length
        for line_idx, line in enumerate( bed ):
            if int( line[ 2 ] ) - int( line[ 1 ] ) == 0:
                print( 'Null length. ' + '\t'.join( line ) + '\nSkip it!' )
                #bed[ line_idx ] = None
                pass
        print( len( bed ) )
        bed = [ line for line in bed if int( line[ 2 ] ) - int( line[ 1 ] ) != 0 ]
        print( len( bed ) )

        ## check in bw with ref in 'chr' format
        if '--add-chr' in sys.argv:
            bed = check_ref_chr( bed, add_chr[ 0 ] )
            pass

        ## retrieve values
        all_inter = retrieve_multi_interval( bed, bw, nbins, group_concat=group_concat )
        print( '\n'.join( [ '\t'.join( [ str( xxx ) for xxx in val_inter ] ) for val_inter in all_inter ] ), file=out_file )
        pass

        if not '---stdout' in sys.argv:
            out_file.close()
            pass
    pass



####
