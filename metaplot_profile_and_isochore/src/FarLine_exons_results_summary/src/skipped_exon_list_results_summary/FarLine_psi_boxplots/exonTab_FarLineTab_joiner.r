#!/usr/bin/Rscript

exonTab_FarLineTab_joiner <- function ( exons_table, FarLine_all_data_tab ) {
  exons_table[ 'join_id' ] <- paste( exons_table$gene_symbol, exons_table$exon_pos, exons_table$exons_flanquants )
  FarLine_all_data_tab[ 'join_id' ] <- paste( FarLine_all_data_tab$gene_symbol, FarLine_all_data_tab$exon_skipped, FarLine_all_data_tab$exons_flanquants )

  merge( exons_table, FarLine_all_data_tab, by='join_id', suffixes=c( '', '.FarLine' ) )
}
