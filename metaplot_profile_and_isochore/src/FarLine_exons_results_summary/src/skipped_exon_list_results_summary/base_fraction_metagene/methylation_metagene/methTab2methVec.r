#!/usr/bin/Rscript

## transform methylation tab in a vector of methylation rate per position
methTab2methVec <- function( meth_tab, chrom, start, end, min.reads=1 ) {
  meth_vec <- rep( NA, end - start + 1 )
  names( meth_vec ) <- as.character( start:end )

  lines_logi <- meth_tab$chrom == chrom & meth_tab$pos >= start & meth_tab$pos <= end
  meth_tab_sub <- meth_tab[ lines_logi, ]

  if ( nrow( meth_tab_sub ) ) {
    for ( row in 1:nrow( meth_tab_sub ) ) {
      xxx <- meth_tab_sub[ row, ]
      if ( xxx$meth_C + xxx$raw_C >= min.reads ) {
        meth_vec[ as.character( xxx$pos ) ] <- xxx$meth_C / ( xxx$meth_C + xxx$raw_C )
      }
    }
  }

  return( meth_vec )
}
