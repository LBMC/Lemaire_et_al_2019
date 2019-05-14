#!/usr/bin/python3

from functions import *
import sys

### MAIN ###
wgbs_bam = sys.argv[ 1 ] #'/home/sebastien/GSE102395_BS_bam_files/MCF7_empytVector_doxWithdrawn_rep1_WGBS_B1s_BestBam_chrSorted.bam'

## load the splicing site sequence table
seq_file = sys.argv[ 2 ] #'/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/FarLine_exons_results_summary/ref_exon_exprim/siDNMT3b-siGL2/results/SS_sequences/inEx50_outOff100/win1b/exon_down_siDNMT3b-siGL2_5SS_ext.seq'
seq_headers, seq_tab = load_file_whead( seq_file )


seq_meth_dict = {}
for seq_id in seq_tab: #range( 1, 5 ):
    ## call the mpileup with samtools, and parse it
    mpileup_dict = mpileup_fc( seq_tab[ seq_id ], wgbs_bam )

    ## recover the position of the CpGs
    CpG_list = CpG_sites_recover( seq_tab[ seq_id ] )

    ## compute the methylation rate
    seq_meth_dict.update( CpG_meth_recover( CpG_list, mpileup_dict ) )
    pass

## write the table of raw/meth cytosines
CpG_meth_out_write( seq_meth_dict )
