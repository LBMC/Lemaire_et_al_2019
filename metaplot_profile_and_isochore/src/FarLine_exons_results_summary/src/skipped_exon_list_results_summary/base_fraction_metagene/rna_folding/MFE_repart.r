#!/usr/bin/Rscript

require( ggplot2 )
require( reshape2 )

all_args <- commandArgs(trailingOnly = F)
MFE_repart_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

## Some function
source( paste( MFE_repart_dir, 'mfe_plots.r', sep='/' ) )
source( paste( MFE_repart_dir, 'rnafold2MFEtable.r', sep='/' ) )
source( paste( MFE_repart_dir, 'mfePerSet.r', sep='/' ) )
source( paste( MFE_repart_dir, 'mfe_stat_test.r', sep='/' ) )


## treat the arguments
args <- commandArgs(trailingOnly = TRUE)
RNAfold_res_vec <- strsplit( args[ 1 ], split=',' )[[ 1 ]]
names_vec <- strsplit( args[ 2 ], split=',' )[[ 1 ]]
out_dir <- args[ 3 ]
dir.create( out_dir, recursive=TRUE, showWarnings=FALSE )
fig_prefix <- paste0( out_dir, '/', args[ 4 ] )

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}


## Recover MFE measure from RNAFold outputs
mfe_df_list <- mfePerSet( RNAfold_res_vec, names_vec=names_vec )

## build the repartition plot
mfe_plots( mfe_df_list, names_vec, fig_prefix, color_pallette=color_pallette )
mfe_stat_test( mfe_df_list, names_vec, stat_test_fc=wilcox.test )
