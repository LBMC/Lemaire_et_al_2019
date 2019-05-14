#!/usr/bin/Rscript

###########
base_counts_fc <- function ( base_tab, ent_base_vec, window_size=1 ) {
  # counts the occurences
  alphabet <- unique( unlist( strsplit( ent_base_vec, split='' ) ) )
  counts <- list()
  for ( letter in ent_base_vec ) {
    letter_vec <- unlist( strsplit( letter, split='' ) )
    counts[[ letter ]] <- apply( base_tab, 2, function( xxx ) {
      sum( xxx %in% letter_vec )
     } )
  }
  counts <- do.call( rbind, counts )

  ## sum counts in the window
  if ( window_size > 1 ) {
    for ( pos in c( 1:( ncol( counts ) - window_size + 1 ) ) ) {
      counts[ , pos ] <- rowSums( counts[ , pos:( pos + window_size - 1 ) ] )
    }
    counts <- counts[ , 1:( ncol( counts ) - window_size + 1 ) ]
  }

  return( counts )
}
