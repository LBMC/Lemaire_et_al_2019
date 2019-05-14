#!/usr/bin/Rscript

fraction_mat_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( fraction_mat_fc_dir, 'fraction_vec_fc.r', sep='/' ) )

###########
fraction_mat_fc <- function ( base_tab, base_vec, window_size, all.pos=TRUE, skew_measure=NULL, use.regexpr=FALSE, presence.only=FALSE ) {
  # adjust right_off according to the lengths of the patterns
  max_base_vec_len <- 1
  if ( is.null( use.regexpr ) ) {
    max_base_vec_len <- max( unlist( lapply( strsplit( base_vec, split='' ), length ) ) )
  }

  ## initiate the matrix of fractions
  fraction_mat <- matrix( -1, nrow=nrow( base_tab ), ncol=ncol( base_tab ) - ( window_size - 1 ) - ( max_base_vec_len - 1 ) )
  left_off <- floor( window_size / 2 )
  right_off <- ceiling( window_size / 2 ) - 1 + ( max_base_vec_len - 1 )

  # adjust left_off and right_off for the setting of the column names in the case of window_size == 1
  if ( left_off == 0 ) {
    left_off <- -1 * ncol( base_tab )
  }
  if ( right_off == 0 ) {
    right_off <- -1 * ncol( base_tab )
  }

  # set the column names
  colnames( fraction_mat ) <- head( tail( colnames( base_tab ), n=-left_off ), n=-right_off )
  rownames( fraction_mat ) <- rownames( base_tab )


  ## compute the fractions
  for ( xxx in 1:dim( base_tab )[ 1 ] ) {
    fraction_mat[ xxx, ] <- fraction_vec_fc( base_tab[ xxx, ], base_vec=base_vec, window_size=window_size, all.pos=all.pos, skew_measure=skew_measure, use.regexpr=use.regexpr, presence.only=presence.only )
  }

  return( fraction_mat )
}
