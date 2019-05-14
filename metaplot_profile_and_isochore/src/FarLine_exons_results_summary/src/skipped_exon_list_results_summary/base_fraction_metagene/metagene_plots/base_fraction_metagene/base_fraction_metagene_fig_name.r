#!/usr/bin/Rscript

require( 'ggplot2', quietly = TRUE )
require( 'reshape2', quietly = TRUE )

## build the curves of metagenes
base_fraction_metagene_fig_name_fc <- function( base_vec, skew_measure=NULL, prefix='', entropy=FALSE, struct_info=FALSE, ent_base_vec=NULL, regexpr_fig_name=NULL, presence.only=FALSE ) {
  # set the name of the figure
  if ( entropy ) {
    if ( is.null( regexpr_fig_name ) ) {
      fig_name <- paste( prefix, paste( ent_base_vec, collapse='-' ), sep='_' )
    } else {
      fig_name <- regexpr_fig_name
    }

    if ( struct_info ) {
      fig_name <- paste( fig_name, 'structInfo.png', sep='_' )
    } else {
      fig_name <- paste( fig_name, 'entropy.png', sep='_' )
    }
  } else {

    if ( is.null( regexpr_fig_name ) ) {
      shaped_base_vec <- unlist( lapply( strsplit( base_vec, split='' ), function( xxx ) { paste( xxx, collapse='p' ) } ) )
      fig_name <- paste( shaped_base_vec, collapse='' )
    } else {
      fig_name <- regexpr_fig_name
    }

    if ( ! is.null( skew_measure ) ) {
      if ( is.null( regexpr_fig_name ) ) {
        fig_name <- paste( fig_name, paste( skew_measure, collapse='' ), sep='-' )
      } else {
        fig_name <- paste0( fig_name, '-skew' )
      }
    }

    if ( presence.only ) {
      fig_name <- paste0( fig_name, '_onlyPres' )
    }

    fig_name <- paste( prefix, fig_name, sep='_' )
    fig_name <- paste0( fig_name, '.png' )
  }

  # message( fig_name )
  return( fig_name )

}
