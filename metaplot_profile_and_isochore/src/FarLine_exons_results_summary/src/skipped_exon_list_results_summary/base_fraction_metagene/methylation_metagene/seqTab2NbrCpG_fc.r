#!/usr/bin/Rscript
require( purrr, quietly=TRUE )

seqTab2NbrCpG_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( seqTab2NbrCpG_fc_dir, 'nb_CpG_fc.r', sep='/' ) )

seqTab2NbrCpG_fc <- function( seq_tab, len.norm=FALSE ) {
  list_seq_vec <- seq_tab$sequence %>% strsplit( split='' )
  lapply( list_seq_vec, nb_CpG_fc, len.norm=len.norm ) %>% return()
}
