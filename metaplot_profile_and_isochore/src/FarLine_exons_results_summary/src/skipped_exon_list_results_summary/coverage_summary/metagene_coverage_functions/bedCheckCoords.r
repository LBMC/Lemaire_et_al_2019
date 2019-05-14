#!/usr/bin/Rscript

#### loader of bed file
bedCheckCoords <- function( annot, clean_dup=FALSE, stranded=FALSE ) {
  annot[ , 1 ] <- as.character( annot[ , 1 ] )
  annot[ , 2 ] <- pmin( annot[[ 2 ]], annot[[ 3 ]] )
  annot[ , 3 ] <- pmax( annot[[ 2 ]], annot[[ 3 ]] )

  annot_clean <- annot
  annot_dup <- NULL
  if ( clean_dup ) {
    if ( stranded ) {
      annot_clean <- annot[ which( ! duplicated( annot[ , c( 1, 2, 3, 6 ) ] ) ), ]
      annot_dup <- annot[ which( duplicated( annot[ , c( 1, 2, 3, 6 ) ] ) ), ]
    } else {
      annot_clean <- annot[ which( ! duplicated( annot[ , c( 1, 2, 3 ) ] ) ), ]
      annot_dup <- annot[ which( duplicated( annot[ , c( 1, 2, 3 ) ] ) ), ]
    }
  }

  return( list( clean=annot_clean, duplicates=annot_dup ) )
}
