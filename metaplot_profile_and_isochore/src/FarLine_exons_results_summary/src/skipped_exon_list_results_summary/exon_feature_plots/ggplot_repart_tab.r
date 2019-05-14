#!/usr/bin/Rscript

ggplot_repart_tab_dir <- dirname(sys.frame(1)$ofile)
source( paste( ggplot_repart_tab_dir, 'repart_fc.r', sep='/' ) )

ggplot_repart_tab <- function( final_data_list, col=NULL, nom_vec=c( 'query' ) ) {

  if ( is.null( col ) ) {
    col <- colnames( final_data_list[[ 1 ]] )[ 1 ]
  }

  list_data_tab <- lapply( final_data_list, function( xxx ) { repart_fc( xxx[ , col ] ) } )

  attr(list_data_tab, "varname") <- "liste"
  data_tab <- melt( list_data_tab, varnames=c( 'value' ), value.name='repart' )
  data_tab$liste <- factor( data_tab$liste, levels=nom_vec )

  return( data_tab )

}
