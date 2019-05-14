#!/usr/bin/Rscript

#### extract one measure from the list of summary tables
measure_extractor <- function( list_total_mm, measure, cond_levels=NULL, diff=FALSE ) {
#   total_mm <- do.call( rbind, lapply( list_total_mm, function(x){
#                          new <- x[ which( x$measure == measure ), ]
#                          return( new )
#   } ) )
  total_mm <- lapply( list_total_mm, function(x){
     new <- x[ which( x$measure == measure ), ]
     return( new )
  } )
  total_mm <- data.frame( do.call( rbind, total_mm ) )

  if ( ! diff ) {
    # condition <- rep( cond_levels, each = dim( total_mm )[ 1 ] / 2 )
    cond_names <- names( list_total_mm )
    nrow_item <- unlist( lapply( total_mm, function (x) { dim( x )[ 1 ] } ) )
    condition <- as.factor( rep( cond_names, nrow_item, levels=cond_levels ) )
    # total_mm <- data.frame( total_mm, condition=rep( condition, each = dim( total_mm )[ 1 ] / length( condition ) ) )
    total_mm <- data.frame( total_mm, condition=condition )
  }
  return( total_mm )
}
