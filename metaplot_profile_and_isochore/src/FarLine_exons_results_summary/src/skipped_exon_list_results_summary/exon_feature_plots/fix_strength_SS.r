#!/usr/bin/Rscript

fix_strength_SS_dir <- dirname(sys.frame(1)$ofile)

fix_strength_SS <- function( feat_tab ) {
  for ( col in c( 'STRENGTH_ACCEPTOR', 'STRENGTH_DONOR', 'FORCE_ACCEPTOR_BEFORE', 'FORCE_ACCEPTOR_AFTER', 'FORCE_DONOR_BEFORE', 'FORCE_DONOR_AFTER') ) {
    idx_bad <- which( feat_tab[ , col ] < 0 )
    feat_tab[ idx_bad, col ] <- NA
  }

  return( feat_tab )
}
