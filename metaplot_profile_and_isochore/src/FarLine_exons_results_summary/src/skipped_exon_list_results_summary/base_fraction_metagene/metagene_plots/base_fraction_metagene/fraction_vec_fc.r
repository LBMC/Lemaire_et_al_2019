#!/usr/bin/Rscript

fraction_vec_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( fraction_vec_fc_dir, 'pattern_finder_fc.r', sep='/' ) )

###########
fraction_vec_fc <- function ( sequence, base_vec, window_size, all.pos=TRUE, skew_measure=NULL, use.regexpr=FALSE, presence.only=FALSE ) {
  # find the bases corresponding to the patterns
  logic_seq <- pattern_finder_fc( sequence, base_vec, use.regexpr=use.regexpr, all.pos=all.pos )
  # print( length( logic_seq ) )

  # compute the fractions, or skewness if not null
  eff_len <- length( logic_seq ) - ( window_size - 1 ) #~ get the length of the output lines
  if ( is.null( skew_measure ) ) {
    fraction_vec <- rep( 0, eff_len )
    for ( pos in c( 1:window_size ) ) {
      fraction_vec <- fraction_vec + logic_seq[ pos:( eff_len + ( pos - 1 ) ) ]
    }
    fraction_vec <- fraction_vec / window_size


    if ( presence.only ) {
      fraction_vec <- as.numeric( as.logical( fraction_vec ) )
    }

  } else { #~ if mode skewness of computation
    counter_logic_seq <- pattern_finder_fc( sequence, skew_measure, use.regexpr=use.regexpr, all.pos=all.pos )

    fraction_vec <- rep( 0, eff_len )
    counter_fraction_vec <- fraction_vec
    for ( pos in c( 1:window_size ) ) {
      fraction_vec <- fraction_vec + logic_seq[ pos:( eff_len + ( pos - 1 ) ) ]
      counter_fraction_vec <- counter_fraction_vec + counter_logic_seq[ pos:( eff_len + ( pos - 1 ) ) ]
    }
    NA_pos <- fraction_vec + counter_fraction_vec == 0
    fraction_vec <- ( fraction_vec - counter_fraction_vec ) / ( fraction_vec + counter_fraction_vec )
    fraction_vec[ NA_pos ] <- NA

  }

  return( fraction_vec )
}
