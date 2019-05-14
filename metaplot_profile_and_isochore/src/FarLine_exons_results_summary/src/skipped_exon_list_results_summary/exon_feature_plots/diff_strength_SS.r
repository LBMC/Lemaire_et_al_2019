#!/usr/bin/Rscript

diff_strength_SS_dir <- dirname(sys.frame(1)$ofile)

diff_strength_SS <- function( feat_tab ) {
  feat_tab[ 'STRENGTH_ACCEPTOR_XvsBEFORE' ] <- feat_tab$STRENGTH_ACCEPTOR - feat_tab$FORCE_ACCEPTOR_BEFORE
  feat_tab[ 'STRENGTH_ACCEPTOR_XvsAFTER' ] <- feat_tab$STRENGTH_ACCEPTOR - feat_tab$FORCE_ACCEPTOR_AFTER
  feat_tab[ 'STRENGTH_DONOR_XvsBEFORE' ] <- feat_tab$STRENGTH_DONOR - feat_tab$FORCE_DONOR_BEFORE
  feat_tab[ 'STRENGTH_DONOR_XvsAFTER' ] <- feat_tab$STRENGTH_DONOR - feat_tab$FORCE_DONOR_AFTER

  return( feat_tab )
}
