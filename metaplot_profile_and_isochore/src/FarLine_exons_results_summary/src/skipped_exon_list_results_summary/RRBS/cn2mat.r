#!/usr/bin/Rscript

cn2mat <- function( cn_tab ) {
  nbr_row <- length( levels( cn_tab$exon_id ) )
  nbr_col <- length( levels( cn_tab$distToSS ) )
  meth_matrix <- rep( NA, nbr_row*nbr_col )

  bbb <- as.numeric( cn_tab$exon_id ) * as.numeric( cn_tab$distToSS )
  meth_matrix[ bbb ] <- cn_tab$methylation_percentage
  meth_matrix <- matrix( meth_matrix, nrow=nbr_row, ncol=nbr_col, dimnames=list( levels( cn_tab$exon_id ), levels( cn_tab$distToSS ) ) )
  return( meth_matrix )
}
