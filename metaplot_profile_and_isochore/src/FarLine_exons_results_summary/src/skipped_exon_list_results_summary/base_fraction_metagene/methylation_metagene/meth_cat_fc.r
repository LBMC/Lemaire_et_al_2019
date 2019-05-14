#!/usr/bin/Rscript

meth_cat_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( meth_cat_fc_dir, 'methTab2methMat.r', sep='/' ) )
source( paste( meth_cat_fc_dir, 'meth_class_fc.r', sep='/' ) )

meth_cat_fc <- function ( meth_tab, seq_tab, min.reads=1 ) {
  ## build matrix of mCpG, base per base
  meth_list <- methTab2methMat( meth_tab, seq_tab, min.reads=min_reads, out_type='list' )
  meth_class <- lapply( meth_list, meth_class_fc, na.rm=T )

  return( meth_class )
}
