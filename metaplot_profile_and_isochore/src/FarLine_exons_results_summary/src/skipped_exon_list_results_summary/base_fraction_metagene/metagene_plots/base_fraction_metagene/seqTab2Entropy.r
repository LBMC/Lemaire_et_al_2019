#!/usr/bin/Rscript

seqTab2Entropy_dir <- dirname(sys.frame(1)$ofile)
source( paste( seqTab2Entropy_dir, 'seq_mat_builder.r', sep='/' ) )
source( paste( seqTab2Entropy_dir, 'base_counts_fc.r', sep='/' ) )
source( paste( seqTab2Entropy_dir, 'entropy_profile_fc.r', sep='/' ) )

###########
seqTab2Entropy <- function ( seq_tab, ent_base_vec, start_pos=1, window_size=1, struct_info=FALSE ) {
  ## build the matrix of sequences
  base_tab <- seq_mat_builder( seq_tab, start_pos=start_pos )

  ## count occurences of each base at each position
  counts_mat <- base_counts_fc( base_tab, ent_base_vec=ent_base_vec, window_size=window_size )

  ## compute the entropy
  entropy_profile <- entropy_profile_fc( counts_mat, struct_info=struct_info )

  return( entropy_profile )

}
