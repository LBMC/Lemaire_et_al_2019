#!/usr/bin/Rscript

ggplot_data_tab <- function ( list_psi_tab, col=NULL, nom_vec=c( 'query' ) ) {

  if ( is.null( col ) ) {
    col <- colnames( list_psi_tab[[ 1 ]] )
  }

  list_data_tab <- lapply( list_psi_tab, function( xxx ) { xxx[ , col ] } )

  data_tab <- melt( list_data_tab )
  data_tab$L1 <- factor( data_tab$L1, levels=nom_vec )
  colnames( data_tab ) <- c( 'variable', 'value', 'liste' )

  return( data_tab )
}
