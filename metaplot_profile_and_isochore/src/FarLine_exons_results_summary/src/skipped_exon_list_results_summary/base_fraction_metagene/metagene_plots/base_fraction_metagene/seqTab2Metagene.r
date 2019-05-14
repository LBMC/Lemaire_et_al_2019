#!/usr/bin/Rscript

seqTab2Metagene_dir <- dirname(sys.frame(1)$ofile)
source( paste( seqTab2Metagene_dir, 'seqTab2FractMat.r', sep='/' ) )
source( paste( seqTab2Metagene_dir, 'mat2Metagene.r', sep='/' ) )

###########
seqTab2Metagene <- function ( seq_tab, start_pos=1, base_vec, reduc_fc=function( xxx ) { mean( xxx, na.rm=TRUE  ) }, window_size, skew_measure=NULL ) {
  ## compute the fractions, proportions
  fraction_mat <- seqTab2FractMat( seq_tab, start_pos=start_pos, base_vec=base_vec, window_size=window_size, skew_measure=skew_measure )

  ## compute the metagene
  metagene <- mat2Metagene( fraction_mat, reduc_fc=reduc_fc )
  return( metagene )
}
