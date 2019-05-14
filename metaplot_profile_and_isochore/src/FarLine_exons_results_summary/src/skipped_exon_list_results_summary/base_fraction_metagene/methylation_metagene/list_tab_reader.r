#!/usr/bin/Rscript

###########
list_tab_reader <- function ( vec_annot_file, vec_annot_name ) {
  list_tab <- list()
  for ( xxx in c( 1:length( vec_annot_file ) ) ) {
    annot_name <- vec_annot_name[ xxx ]
    annot_file <- vec_annot_file[ xxx ]
    list_tab[[ annot_name ]] <- read.table( annot_file, stringsAsFactors=FALSE, sep='', header=TRUE )
  }
  return( list_tab )
}
