#!/usr/bin/Rscript

###########
entropy_fc <- function ( counts_vec ) {
  fractions <- counts_vec / sum( counts_vec )
  terms <- fractions * log( fractions )
  terms[ counts_vec == 0 ] <- 0
  entropy_val <- -sum( terms )

  return( entropy_val )

}
