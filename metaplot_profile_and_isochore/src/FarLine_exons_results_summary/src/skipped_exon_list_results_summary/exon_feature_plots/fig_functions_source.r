#!/usr/bin/Rscript

fig_functions_source <- function( parent_path ) {
  source( paste( parent_path, 'ggplot_data_tab.r', sep='/' ) )
  source( paste( parent_path, 'ggplot_repart_tab.r', sep='/' ) )
  source( paste( parent_path, 'fig_width.r', sep='/' ) )
  source( paste( parent_path, 'fig_name.r', sep='/' ) )
}
