#!/usr/bin/python3

import sys

def dup_rem( liste ):
    items = set()
    for idx, xxx in enumerate( liste ):
        if xxx in items:
            liste[ idx ] = None
        else:
            items.add( xxx )
            pass

    liste = [ xxx for xxx in liste if xxx != None ]
    return( liste )

ref_order_path = sys.argv[ 1 ]
# ref_order_path = './report_cov_lists/histoneMeth_complete_refined'
ratio_tab = sys.argv[ 2 ]
# ratio_tab = './summary_mean_cov_ratio.tsv'

with open( ref_order_path, 'r' ) as rofile:
    ref_order_list = rofile.readlines()
    ref_order_list = [ xxx.strip() for xxx in ref_order_list if xxx[ 0 ] not in [ '#', '\n' ] ]
    pass
ref_order_list = dup_rem( ref_order_list )


with open( ratio_tab, 'r' ) as ratio_file:
    ratio_h = ratio_file.readline().strip().split( '\t' )
    ed_idx = ratio_h.index( 'exp_des' )

    stats_dic = {}
    for idx, line in enumerate( ratio_file ):
        parts = line.strip().split( '\t' )
        if parts[ ed_idx ] in stats_dic:
            stats_dic[ parts[ ed_idx ] ].append( line )
        else:
            stats_dic[ parts[ ed_idx ] ] = [ line ]
            pass
        pass
    pass


print( '\t'.join( ratio_h ), file=sys.stdout )
for xxx in ref_order_list:
    if xxx in stats_dic:
        for line in stats_dic[ xxx ]:
            print( line, end='', file=sys.stdout )
            pass
        pass
    else:
        print( 'No ratio for: %s' % xxx, file=sys.stderr )
        pass
    pass


####
