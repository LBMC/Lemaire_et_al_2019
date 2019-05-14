#!/usr/bin/Rscript

meth_meas_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( meth_meas_fc_dir, 'window_summary.r', sep='/' ) )
source( paste( meth_meas_fc_dir, 'methTab2methMat.r', sep='/' ) )

meth_meas_fc <- function ( meth_tab, seq_tab, min.reads=1, start_pos=1, window=1 ) {
  ## build matrix of mCpG, base per base
  meth_list <- methTab2methMat( meth_tab, seq_tab, min.reads=min_reads, start_pos=start_pos, out_type='list' )
  meth_meas <- lapply( meth_list, mean, na.rm=T )

  return( meth_meas )
}
