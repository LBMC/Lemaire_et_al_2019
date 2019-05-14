#!/usr/bin/Rscript

## build the curves of metagenes
sub_col_pal_fc <- function( query_mat, scale_range=c( 0, 1 ), color_pallette=heat.colors( 256 ) ) {
  # range of values for query matrix
  query_range <- range( query_mat, na.rm=TRUE )

  # relative query range to the scale range
  rel_query_range <- query_range
  for ( xxx in 1:length( query_range ) ) {
    rel_query_range[ xxx ] <- ( query_range[ xxx ] - scale_range[ 1 ] ) / ( scale_range[ 2 ] - scale_range[ 1 ] )
  }

  # border idx in color_pallette for query values
  col_pal_idx_border <- round( quantile( c( 1:length( color_pallette ) ), probs=rel_query_range ) )

  # color_pallette for query values
  query_col_pal <- color_pallette[ col_pal_idx_border[ 1 ]:col_pal_idx_border[ 2 ] ]

  return( query_col_pal )
}
