#!/usr/bin/Rscript

###########
seq_mat_builder <- function ( seq_tab, start_pos=1 ) {
  base_tab <- do.call( rbind, strsplit( seq_tab$sequence, split='' ) )
  colnames( base_tab ) <- ( start_pos:( start_pos + ncol( base_tab ) -1 ) )
  rownames( base_tab ) <- paste( seq_tab$coordinates, seq_tab$strand, sep=':' )
  return( base_tab )
}
