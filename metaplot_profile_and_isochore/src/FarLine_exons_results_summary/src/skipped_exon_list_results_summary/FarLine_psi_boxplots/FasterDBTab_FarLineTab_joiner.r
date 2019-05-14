#!/usr/bin/Rscript

FasterDBTab_FarLineTab_joiner <- function ( FasterDB_feature_tab, FarLine_all_data_tab ) {
  FasterDB_feature_tab[ 'join_id' ] <- paste( FasterDB_feature_tab$GENE, FasterDB_feature_tab$EXON )
  FarLine_all_data_tab[ 'join_id' ] <- paste( FarLine_all_data_tab$id_gene, FarLine_all_data_tab$exon_skipped )

  merge( FasterDB_feature_tab, FarLine_all_data_tab, by='join_id', suffixes=c( '.FDB', '.FarLine' ) )
}
