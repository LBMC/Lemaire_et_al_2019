#!/usr/bin/Rscript

repart_fc <- function( xxx ) {
  occ <- table( xxx )
  cum_occ <- cumsum( occ )
  repart <- cum_occ / tail( cum_occ, n=1 )
  repart <- as.table( repart )

  return( repart )
}
