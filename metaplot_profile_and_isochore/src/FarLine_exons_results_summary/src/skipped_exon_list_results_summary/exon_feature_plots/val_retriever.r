#!/usr/bin/Rscript

# retrieve value from gathered data in one string
val_retriever <- function ( data_vec, sep='', field=1, stringsAsFactors=FALSE, data_type=as.character ) {
  data_type( as.data.frame( do.call( rbind, strsplit( data_vec, split=sep ) ), stringsAsFactors=stringsAsFactors )[ , field ] )
}
