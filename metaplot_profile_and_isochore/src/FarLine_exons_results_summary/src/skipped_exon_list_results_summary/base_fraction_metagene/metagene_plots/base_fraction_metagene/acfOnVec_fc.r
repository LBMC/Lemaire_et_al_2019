#!/usr/bin/Rscript

acfOnVec_fc <- function( vec ) {
  vec_acf <- acf( vec, plot=FALSE, na.action=na.pass )[[ 'acf' ]][ , 1, 1 ]

  return( vec_acf )
}
