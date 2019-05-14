#!/usr/bin/Rscript

nb_CpG_fc_dir <- dirname(sys.frame(1)$ofile)

nb_CpG_fc <- function( seq, len.norm=FALSE ) {
  # seq <- seq_str %>% strsplit( split='' ) %>% unlist()
  nbr_CpG <- paste0( head( seq, n=-1 ), tail( seq, n=-1 ) ) %in% 'CG' %>% sum()
  if ( len.norm ) {
    nbr_CpG <- nbr_CpG / ( length( seq ) - 1 )
  }

  return( nbr_CpG )
}
