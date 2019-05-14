#!/usr/bin/Rscript

###########
regexpr_finder_fc <- function ( sequence, pattern, all.pos=TRUE ) {
    seq_str <- paste( sequence, collapse='' )

    # find the matches
    matches <- gregexpr(
            pattern=pattern,
            text=seq_str,
            perl=TRUE,
          )[[ 1 ]]


    # record the base positions of the matches
    log_vec <- rep( FALSE, nchar( seq_str ) )

    if ( matches[ 1 ] != -1 ) {
      for ( xxx in c( 1:length( matches ) ) ) {
        start <- matches[ xxx ]

        end <- start
        if ( all.pos ) {
          end <- start + attr( matches, 'match.length' )[ xxx ] - 1
        }
        log_vec[ start:end ] <- TRUE
      }
    }


  return( log_vec )
}
