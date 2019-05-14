#!/usr/bin/Rscript

seqTab2FractMat_dir <- dirname(sys.frame(1)$ofile)
source( paste( seqTab2FractMat_dir, 'seq_mat_builder.r', sep='/' ) )
source( paste( seqTab2FractMat_dir, 'fraction_mat_fc.r', sep='/' ) )

###########
seqTab2FractMat <- function ( seq_tab, start_pos=1, base_vec, window_size, all.pos=TRUE, skew_measure=NULL, use.regexpr=FALSE, presence.only=FALSE ) {
  ## build the matrix of sequences
  base_tab <- seq_mat_builder( seq_tab, start_pos=start_pos )

  ## compute the fractions, proportions
  fraction_mat <- fraction_mat_fc( base_tab, base_vec=base_vec, window_size=window_size, all.pos=all.pos, skew_measure=skew_measure, use.regexpr=use.regexpr, presence.only=presence.only )
  return( fraction_mat )
}
