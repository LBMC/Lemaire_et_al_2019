#!/usr/bin/Rscript

fig_name <- function( plot_name, out_dir, nom, fig_type ) {
  paste( out_dir, '/', nom, '_', plot_name, '_', fig_type, '.png', sep='' )
}
