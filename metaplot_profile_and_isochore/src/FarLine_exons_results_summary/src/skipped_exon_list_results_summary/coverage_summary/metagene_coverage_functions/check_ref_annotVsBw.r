#!/usr/bin/Rscript

#### check refs between BigWig file and annotations
check_ref_annotVsBw <- function( annot, bw_file, ref_add_chr=F ) {
  bw_reader <- load.bigWig( bw_file )

  annot_ref_col <- annot[ , 1 ]
  if ( ref_add_chr ) {
    annot_ref_col <- paste0( 'chr', annot_ref_col )
  }
  annot_ref <- unique( annot_ref_col )

  mis_ref <- setdiff( annot_ref, bw_reader$chroms )
  if ( length( mis_ref ) == length( annot_ref ) ) {
    message( 'All ref missing!' )
    message( paste( annot[ 1, ], collapse='\t' ) )
    stop()
  } else if ( length( mis_ref ) != 0 ) {
    message( 'Several ref missing!' )
    message( paste( mis_ref, collapse='\t' ) )
    message( 'Filter the annotations!' )
    annot <- annot[ which( ! sapply( annot_ref_col, function( xxx, yyy ) { xxx %in% yyy }, mis_ref ) ), ]
  }

  return( annot )
}
