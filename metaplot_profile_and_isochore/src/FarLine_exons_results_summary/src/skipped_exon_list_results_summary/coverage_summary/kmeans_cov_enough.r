#!/usr/bin/Rscript

kmeans_cov_enough <- function ( comp_pair, cov_df, min_max_signal ) {
  cov_enough_id <- c( 1:nrow( cov_df ) )

  if ( length( comp_pair ) == 1 ) {
    cov_enough_id <- intersect( cov_enough_id, which( apply( cov_df, 1, max ) > min_max_signal ) )
  } else if ( length( comp_pair ) == 2 ) {
    cov_enough_id <- intersect( cov_enough_id, which( apply( abs( cov_df ), 1, max ) > -1 ) )
  }
}
