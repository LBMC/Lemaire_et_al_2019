#!/usr/bin/python3

import sys, re

def out_tab_line( list_item, channel=sys.stdout ):
    print( '\t'.join( list_item ), file=channel )
    pass

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

def ratio_fc( nume, denom ):
    if denom:
        ratio_val = nume / denom
    else:
        ratio_val = 'NA'
        pass

    return( ratio_val )

def diff_norm_fc( nume, denom ):
    if denom:
        ratio_val = ( nume - denom ) / denom
    else:
        ratio_val = 'NA'
        pass

    return( ratio_val )

def diff_norm_max_fc( nume, denom ):
    if nume or denom:
        ratio_val = ( nume - denom ) / max( nume, denom )
    else:
        ratio_val = 'NA'
        pass

    return( ratio_val )

####

print( sys.argv, file=sys.stderr )

smc_path = sys.argv[ 1 ]
# smc_path = '/home/slemaire/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/res/geneVisu_exon/results/figures/coverage_summary/varLen1000/ENCSR000AKC_ENCFF206QIK_H3K27ac_GM12878_1/cov/summary_mean_cov.tsv'
ratio_list = [ xxx.split( '/' ) for xxx in sys.argv[ 2 ].split( ',' ) ]
# ratio_list_arg = 'GC_rich_exon/CCE_exon,AT_rich_exon/CCE_exon,GC_rich_exon/AT_rich_exon'
# ratio_list = [['GC_rich_exon', 'CCE_exon'], ['AT_rich_exon', 'CCE_exon'], ['GC_rich_exon', 'AT_rich_exon']]

ratio_type_dic = { 'r': ratio_fc, 'dn': diff_norm_fc, 'dnm': diff_norm_max_fc }
ratio_type_list = [ 'r' for xxx in ratio_list ]
flag = '--ratio-type'
if flag in sys.argv:
    ratio_type_list = sys.argv[ sys.argv.index( flag ) + 1 ].split( ',' )
    pass

with open( smc_path ) as smc_file:
    smc_lines = [ xxx[ :-1 ].split( '\t' ) for xxx in smc_file.readlines() ]
    smc_h = smc_lines[ 0 ]
    smc_lines = smc_lines[ 1: ]
    pass


smc_h = smc_h + [ 'sample' ]
cond_idx = smc_h.index( 'cond' )
rep_idx = smc_h.index( 'rep' )
samp_idx = smc_h.index( 'sample' )
if 'exp_des' in smc_h:
    exp_des_idx = smc_h.index( 'exp_des' )
    pass
for line_idx, line in enumerate( smc_lines ):
    exp_des_str = ''
    if 'exp_des' in smc_h:
        exp_des_str = '_' + smc_lines[ line_idx ][ exp_des_idx ]
        pass
    smc_lines[ line_idx ] = smc_lines[ line_idx ] + [ smc_lines[ line_idx ][ cond_idx ] + '_' + smc_lines[ line_idx ][ rep_idx ] + exp_des_str ]
    pass

# print( smc_h, file=sys.stderr )
# print( smc_lines[ 0 ], file=sys.stderr )

samp_list = dup_rem( [ xxx[ samp_idx ] for xxx in smc_lines ] )

annot_idx = smc_h.index( 'annot_list' )
out_tab_line( [ 'measure', 'cond', 'rep', 'exp_des' ] + [ '/'.join( xxx ) for xxx in ratio_list ] )
for measure in smc_h[ 0:6 ]:
    meas_idx = smc_h.index( measure )
    # print( '--- ' + measure, file=sys.stdout )
    for samp in samp_list:
        # print( '--- ' + samp, file=sys.stderr )
        smc_lines_sub = [ xxx for xxx in smc_lines if xxx[ samp_idx ] == samp ]
        # print( smc_lines_sub, file=sys.stderr )
        if samp == 'GSE75792_GSM1967871_K562_POLR2AphosphoS2_B1s':
            print( smc_lines_sub, file=sys.stderr )
            pass

        ratio_val_list = [0.0]*len( ratio_list )
        for ratio_idx, ratio in enumerate( ratio_list ):
            try:
                nume = float( [ xxx[ meas_idx ] for xxx in smc_lines_sub if re.search( pattern=ratio[ 0 ], string=xxx[ annot_idx ] ) ][ 0 ] )
                denom = float( [ xxx[ meas_idx ] for xxx in smc_lines_sub if re.search( pattern=ratio[ 1 ], string=xxx[ annot_idx ] ) ][ 0 ] )
            except:
                if len( smc_lines_sub ) == 0:
                    print( 'Missing ChIP sample: %s' % ( samp ), file=sys.stderr )
                    pass
                else:
                    annot_sub = [ xxx[ annot_idx ] for xxx in smc_lines_sub ]
                    diff_set = set( ratio ).difference( set( annot_sub ) )
                    print( 'Missing partners: %s, for ratio in samp: %s' % ( samp, ' '.join( list( diff_set ) ) ), file=sys.stderr )
                    pass
                ratio_val_list[ ratio_idx ] = 'X'
                continue


            ratio_val_list[ ratio_idx ] = ratio_type_dic[ ratio_type_list[ ratio_idx ] ]( nume, denom )
            pass

        out_tab_line( [ measure, smc_lines_sub[ 0 ][ cond_idx ], smc_lines_sub[ 0 ][ rep_idx ], smc_lines_sub[ 0 ][ exp_des_idx ] ] + [ str( xxx ) for xxx in ratio_val_list ] )
        pass
    pass




####
