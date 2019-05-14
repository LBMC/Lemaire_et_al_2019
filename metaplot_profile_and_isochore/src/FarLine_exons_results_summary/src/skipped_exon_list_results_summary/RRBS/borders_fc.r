#!/usr/bin/Rscript

borders_fc <- function( limits ) {
  borders <- data.frame( basse=c( -Inf, limits ), haute=c( limits, Inf ) )
  rownames( borders ) <- seq( 1, nrow( borders ) )

  # build labels of the categories
  domains <- paste("[",":", borders$haute[1], "]", sep='')
  if (length(limits) > 1) {
    for (xxx in seq(2, length(limits))) {
      domains <- c(domains, paste("]", borders$basse[xxx], ":", borders$haute[xxx], "]", sep=''))
    }
  }
  domains <- c(domains, paste("]", tail(borders$basse, n=1), ":", "]", sep=''))
  domains <- factor( domains, levels=domains )

  borders[ 'inter' ] <- domains
  return( borders )
}
