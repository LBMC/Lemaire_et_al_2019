#!/usr/bin/Rscript

meth_class_fc_dir <- dirname(sys.frame(1)$ofile)

meth_class_fc <- function( meth_line, na.rm=TRUE, fraction=TRUE ) {
  meth_class_list <- list(
    'Meth'=sum( meth_line >= 0.82, na.rm=na.rm ),
    'Inter'=sum( meth_line < 0.82 & meth_line > 0.13, na.rm=na.rm ),
    'Raw'=sum( meth_line <= 0.13, na.rm=na.rm )
  )

  if ( fraction ) {
    nbr_CpG_eval <- sum( ! is.na( meth_line ) )
    meth_class_list <- lapply( meth_class_list, function( xxx ) { xxx / nbr_CpG_eval } )
  }

  return( meth_class_list )
}
