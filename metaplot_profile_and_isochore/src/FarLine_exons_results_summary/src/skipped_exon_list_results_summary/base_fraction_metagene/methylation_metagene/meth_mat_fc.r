#!/usr/bin/Rscript

meth_mat_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( meth_mat_fc_dir, 'window_summary.r', sep='/' ) )
source( paste( meth_mat_fc_dir, 'methTab2methMat.r', sep='/' ) )

meth_mat_fc <- function ( meth_tab, seq_tab, min.reads=1, start_pos=1, window=1 ) {
  ## build matrix of mCpG, base per base
  meth_mat <- methTab2methMat( meth_tab, seq_tab, min.reads=min_reads, start_pos=start_pos )

  ## average the methylation rate in a scanning window
  meth_mat <- window_summary( meth_mat, window=window )

  return( meth_mat )
}
