#!/usr/bin/Rscript

## compute meta-methylation rate per position from a matrix of methylation rate of a set of annotation
meth_meta_fc <- function( meth_mat ) {
  ## metaplot of the methylations
  meth_meta <- colMeans( meth_mat, na.rm=TRUE )
  # meth_meta[ is.nan( meth_meta ) ] <- NA

  ## count nbr of exons with CpG at each position
  CpGperPos <- apply( meth_mat, 2, function( xxx ) { ( length( xxx ) - sum( is.na( xxx ) ) ) / length( xxx ) } )

  return( list( 'meth_meta'=meth_meta, 'CpG_counts'=CpGperPos ) )
}
