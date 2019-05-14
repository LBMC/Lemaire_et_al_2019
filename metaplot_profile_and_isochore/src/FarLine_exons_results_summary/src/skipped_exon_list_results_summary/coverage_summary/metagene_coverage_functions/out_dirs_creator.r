#!/usr/bin/Rscript

#### create the output folders and return their paths
out_dirs_creator <- function ( out_dir, all=TRUE ) {
  list_out_dir <- list()
  list_out_dir[[ 'metagene_mean' ]] <- paste( out_dir, 'metagene_mean', sep='/' );
  list_out_dir[[ 'metagene_median' ]] <- paste( out_dir, 'metagene_median', sep='/' );
  list_out_dir[[ 'heatmapSimple' ]] <- paste( out_dir, 'heatmapSimple', sep='/' );
  list_out_dir[[ 'repartMax' ]] <- paste( out_dir, 'repartMax', sep='/' );
  list_out_dir[[ 'repartMean' ]] <- paste( out_dir, 'repartMean', sep='/' );
  if ( all ) {
    list_out_dir[[ 'metagene_summary' ]] <- paste( out_dir, 'metagene_summary', sep='/' );
    list_out_dir[[ 'heatmapNormalized' ]] <- paste( out_dir, 'heatmapNormalized', sep='/' );
  }

  lapply( list_out_dir, function ( xxx ) { dir.create( file.path( xxx ), showWarnings = FALSE, recursive=TRUE ) } )

  return( list_out_dir )
}
