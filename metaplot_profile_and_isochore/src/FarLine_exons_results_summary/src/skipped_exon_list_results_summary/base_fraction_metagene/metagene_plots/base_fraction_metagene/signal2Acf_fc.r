#!/usr/bin/Rscript

signal2Acf_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( signal2Acf_fc_dir, 'per_mat_fc.r', sep='/' ) )

signal2Acf_fc <- function( fraction_mat ) {
  per_mat <- per_mat_fc( fraction_mat )[[ 'per_mat' ]]
  per_meta <- colMeans( per_mat, na.rm=TRUE )

  return( per_meta )
}
