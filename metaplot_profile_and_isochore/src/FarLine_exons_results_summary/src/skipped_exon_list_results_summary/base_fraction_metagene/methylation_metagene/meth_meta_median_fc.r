#!/usr/bin/Rscript

require( tidyverse, quietly=TRUE )
require( reshape2 )

## compute median among the replicate
meth_meta_median_fc <- function( all_list_meth_meta ) {
  # add column mith rownames ( relative base positions )
  dada <- modify_depth( .x=all_list_meth_meta, .depth=4, .ragged=TRUE, .f=function( xxx ) {
    data.frame( value=xxx, idx=names( xxx ), stringsAsFactors=FALSE )
  } )

  # convert in data.frame
  clac <- melt( dada, id.vars=c( 'idx', 'value' ) )
  clac$idx <- as.numeric( clac$idx )

  # compute the median among the replicates
  rep_lev <- unique( clac$L2 )

  meth_meta_df <- clac[ clac$L4 == "meth_meta", !( names( clac ) == 'L4' ) ]
  meth_meta_df <- spread( meth_meta_df, L2, value )
  meth_meta_df[ 'median' ] <- apply( meth_meta_df[ rep_lev ], 1, median, na.rm=TRUE )

  # recover the CpG_counts corresponding to the median value
  CpG_counts_df <- clac[ clac$L4 == "CpG_counts", !( names( clac ) == 'L4' ) ]
  CpG_counts_df <- spread( CpG_counts_df, L2, value )

  CpG_log_mat <- as.matrix( meth_meta_df[ rep_lev ] ) == meth_meta_df$median
  CpG_log_mat[ is.na( CpG_log_mat ) ] <- FALSE

  CpG_mat <- matrix( NA, nrow=length( meth_meta_df$median ), ncol=length( rep_lev ) )
  CpG_mat[ CpG_log_mat ] <- as.matrix( CpG_counts_df[ rep_lev ] )[ CpG_log_mat ]
  meth_meta_df[ 'CpG_counts' ] <- apply( CpG_mat, 1, function( xxx ) {
    if ( all( is.na( xxx ) ) ) {
      return( 0 )
    } else {
      return( max( xxx, na.rm=TRUE ) )
    }
  })

  # re-convert in list
  meth_meta_df <- meth_meta_df[ , !( names( meth_meta_df ) %in% rep_lev ) ]
  meth_meta_list <- lapply( split( meth_meta_df[ , names( meth_meta_df ) != 'L1' ], meth_meta_df$L1, drop = TRUE ), function( x ) split( x[ , names( x ) != 'L3' ], x[[ 'L3' ]], drop = TRUE ) )

  return( meth_meta_list )
}
