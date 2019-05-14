#!/usr/bin/Rscript

fraction_computer <- function( meth_matrix, borders, start_pos ) {
  # counts the number of cytosine, with methylation in a range, for each position (column)
  counts <- list()
  counts[[ 'tot_eval' ]] <- colSums( !is.na( meth_matrix ) )
  for ( xxx in seq( 1, nrow( borders ) ) ) {
    counts[[ as.character( borders$inter[ xxx ] ) ]] <- colSums( meth_matrix > borders$basse[ xxx ] & meth_matrix <= borders$haute[ xxx ], na.rm=TRUE )
  }
  counts <- do.call( rbind, counts )

  # sum the counts in a sliding window (final count being at the center of the window)
  window_counts <- matrix( 0,
    nrow=nrow( counts ),
    ncol=ncol( counts ) - ( window - 1 ),
    dimnames=list( rownames( counts ),
      colnames( counts )[ 1:( ncol( counts ) - ( window - 1 ) ) ]
      )
    )
  for ( xxx in seq( 1, window ) ) {
    window_counts <- window_counts + counts[ , xxx:( xxx + ( ncol( window_counts ) - 1 ) ) ]
  }
  colnames( window_counts ) <- seq( start_pos, start_pos + ncol( window_counts ) - 1 )

  # compute the fractions of each meth-cytosine category in each window
  null_loc <- window_counts[ 1, ] == 0
  fractions <- t( apply( window_counts[ -1, ], 1, function (x) { x[ !null_loc ] <- x[ !null_loc ] / window_counts[ 1, !null_loc ]; x[ null_loc ] <- NA; return( x ) } ) )

  return( fractions )
}
