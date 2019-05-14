#!/usr/bin/Rscript

#### bed file loader
bed_loader <- function( annot_file, clean_dup=FALSE, stranded=FALSE ) {
  annot <- read.table( annot_file, header=FALSE, stringsAsFactors=FALSE )
  if ( length( grep( 'ppp', annot[ , 5 ] ) ) > 0 ) {
      annot <- cbind( annot, strand=do.call( rbind, strsplit( annot[ ,5 ], split='ppp' ) )[ ,2 ] )
  }
  annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )
  return( annot )
}
