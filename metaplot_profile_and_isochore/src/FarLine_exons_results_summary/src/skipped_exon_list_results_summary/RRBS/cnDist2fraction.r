#!/usr/bin/Rscript

all_args <- commandArgs(trailingOnly = F)
cnDist2fraction_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( cnDist2fraction_dir, 'borders_fc.r', sep='/' ) )
source( paste( cnDist2fraction_dir, 'cn2mat.r', sep='/' ) )
source( paste( cnDist2fraction_dir, 'fraction_computer.r', sep='/' ) )
source( paste( cnDist2fraction_dir, 'reg_ext_window.r', sep='/' ) )

tabLoader <- function ( cn_input, start_pos=1, end_pos=NULL, window_size=1 ) {
  raw <- read.table( cn_input, header=TRUE )

  if ( is.null( end_pos ) ) {
    end_pos=as.numeric( max( raw$pos ) )
  }
  range_pos <- seq( start_pos - floor( window_size / 2 ), end_pos + ceiling( window_size / 2 ) - 1 )

  raw$exon_id <- as.factor( raw$exon_id )
  raw$distToSS <- factor( raw$distToSS , levels=as.character( range_pos ) )
  return( raw )
}

#########

args <- commandArgs(trailingOnly = TRUE)
cn_input_vec <- unlist( strsplit( args[1], split=',' ) ) #~ '~/analyses_SEBASTIEN/analyse_RRBS_auboeuf_MCF7-SH5Y/chrSorted_bismarkBam/CG_context/distToSplicingSites_cytosineMethylationRate/CN_files/CG_RRBS_MCF7siGL2n1_5mCrate_distTo3pSS_exons_up_1kbUpDown3pSS.cn'
limits <- as.numeric( unlist( strsplit( args[2], split=',' ) ) ) #~ c( 13, 87 )
rel_pos <- as.numeric( unlist( strsplit( args[3], split=',' ) ) )
window <- as.numeric( args[4] ) #~ 10
out_file <- args[5]
conditions <- unlist( strsplit( args[6], split=',' ) )
replicates <- unlist( strsplit( args[7], split=',' ) )

# #########
# ## compute range of positions
# range_pos <- as.factor( seq( rel_pos[1], rel_pos[2] ) )

#########
## build the matrix of methylation rate per annotation and per position
df_expdesign <- data.frame( file=cn_input_vec, condition=conditions, replicate=replicates, row.names=seq( 1, length( cn_input_vec ) ) )
raw_list <- lapply( cn_input_vec, tabLoader, start_pos=rel_pos[ 1 ], end_pos=rel_pos[ 2 ], window_size=window )
meth_matrix_list <- lapply( raw_list, cn2mat )

#########
borders <- borders_fc( limits )

#########
fractions_list <- lapply( meth_matrix_list, fraction_computer, borders=borders, start_pos=rel_pos[ 1 ] )

#########
saveRDS( list( fractions=fractions_list, df_expdesign=df_expdesign ), out_file )
