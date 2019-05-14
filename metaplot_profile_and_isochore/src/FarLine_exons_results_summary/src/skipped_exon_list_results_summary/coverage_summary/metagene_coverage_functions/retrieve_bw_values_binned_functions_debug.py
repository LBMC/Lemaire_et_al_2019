#!/usr/bin/python3

## Usage: python3 <script> <bw1,bw2,...> <bed1,bed2,...> <out_dir> <nbins> [ --add-chr 0,1,... ]

import sys, os
import pyBigWig
from math import isnan
from statistics import mean
import operator


####
def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))

def load_table( file_path ):
    file_open = open( file_path, 'r' )
    tab_lines = [ xxx.strip().split( '\t' ) for xxx in file_open.readlines() ]
    file_open.close()
    return( tab_lines )

def get_items_unique( tab, col ):
    group_list = [ '' ]*len( tab )
    for idx, xxx in enumerate( tab ):
        if xxx[ col ] not in group_list:
            group_list[ idx ] = xxx[ col ]
            pass
        pass
    group_list = [ xxx for xxx in group_list if xxx != '' ]

    return( group_list )

def val_inter_NoneZero( ref, start, end, bw ):
    intervalle = bw.values( ref, start, end )
    nan_val = float( 'nan' )
    for idx, xxx in enumerate( intervalle ):
        if xxx == None or isnan( xxx ):
            intervalle[ idx ] = 0
            pass
        pass
    return( intervalle )

def check_ref_chr( bed, add_chr ):
    if add_chr and not bed[ 0 ][ 0 ].startswith( 'chr' ):
        for annot in bed:
            annot[ 0 ] = 'chr' + annot[ 0 ]
            pass
    elif not add_chr and bed[ 0 ][ 0 ].startswith( 'chr' ):
        for annot in bed:
            annot[ 0 ] = annot[ 0 ][ 3: ]
            pass
        pass

    return( bed )

def bin_convert( vec, nbins ):
    new_inter = [None]*nbins
    work_len = len( vec ) - 1
    ends = [ int( ( float( xxx + 1 ) / nbins ) * work_len ) for xxx in range( nbins ) ]
    starts = [ 0 ] + [ xxx for xxx in ends[ :-1 ] ]
    for pos in range( nbins ):
        new_inter[ pos ] = mean( vec[ starts[ pos ]:( ends[ pos ] + 1 ) ] )
        pass
    return( new_inter )

def retrieve_val_inter( annot, bw, nbins=None ):
    ref = annot[ 0 ]
    start = int( annot[ 1 ] )
    end = int( annot[ 2 ] )
    try:
        val_inter = val_inter_NoneZero( ref, start, end, bw )
    except RuntimeError:
        print( 'RuntimeError: Bad interval in val_inter_NoneZero function' )
        print( ( ref, start, end, bw ) )
        sys.exit()

    ## check strand
    if annot[ 5 ] == '-1':
        val_inter.reverse()
        pass

    ## convert in bins
    if nbins != None:
        try:
            val_inter = bin_convert( val_inter, nbins )
        except:
            print( len( val_inter ) )
            print( annot )
        pass

    return( val_inter )

def retrieve_grouped_val_inter( annot_list, bw, nbins=None ):
    val_inter = []
    for annot in annot_list:
        val_inter = val_inter + retrieve_val_inter( annot, bw, nbins=None )
        pass

    ## convert in bins
    if nbins != None:
        try:
            val_inter = bin_convert( val_inter, nbins )
        except:
            print( len( val_inter ) )
            print( annot )
        pass

    return( val_inter )

def retrieve_multi_interval( bed, bw, nbins, group_concat=False ):
    if not group_concat:
        all_inter = [None]*len( bed )
        for annot_idx, annot in enumerate( bed ):
            all_inter[ annot_idx ] = retrieve_val_inter( annot, bw, nbins=nbins )
            # print( '\t'.join( [ str( xxx ) for xxx in all_inter[ annot_idx ] ] ), file=sys.stdout )
            pass
    else:
        bed = [ xxx + xxx[ 3 ].split( '_' ) for xxx in bed ]
        for xxx in bed:
            xxx[ 7 ] = int( xxx[ 7 ] )
            pass

        ## get the groups of annotations
        group_list = get_items_unique( bed, 6 )
        all_inter = [None]*len( group_list )
        for group_idx, group in enumerate( group_list ):
            sub_tab = [ xxx for xxx in bed if xxx[ 6 ] == group ]
            sub_tab = sort_table( sub_tab, 7 )
            all_inter[ group_idx ] = retrieve_grouped_val_inter( sub_tab, bw, nbins=nbins )
            pass

    return( all_inter )


def retrieve_multi_mean_cov( bed, bw ):
    all_means = [None]*len( bed )
    for annot_idx, annot in enumerate( bed ):
        # print( str( annot_idx ) + '/' + str( len( bed ) ), end='\r', file=sys.stderr )
        inter = retrieve_val_inter( annot, bw, nbins=None )
        all_means[ annot_idx ] = sum( inter ) / len( inter )
        pass
    # print( '', file=sys.stderr )

    return( all_means )
####
