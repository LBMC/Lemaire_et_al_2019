#!/usr/bin/Rscript

require( 'bigWig' )

## absolute coverage in one condition
stranded_coverages <- function( bw_file, annot, off_set=1, normFactor=1, harmonization = FALSE, ref_add_chr=F ) {
  bw_reader <- load.bigWig( filename=bw_file )

  ## check if human chromosome
  if ( ref_add_chr ) {
    annot[ , 1 ] <- paste0( 'chr', annot[ , 1 ] )
  }

  # retrieve coverage from bam file for all annotations and apply normalization
  #~ cov_df <- bed.step.probeQuery.bigWig( bw=bw_reader, bed=annot, step=1, as.matrix=TRUE, gap.value=0 ) / normFactor
  cov_df <- bed.step.probeQuery.bigWig( bw=bw_reader, bed=annot, step=1, as.matrix=TRUE, gap.value=0 )
  if ( ! harmonization ) {
    cov_df <- cov_df / normFactor
  }
  #~ cov_df[ is.na( cov_df ) ] <- 0

  cov_df <- data.frame( cov_df )
  if ( dim( annot )[ 2 ] >= 6 ) {
    cov_df[ annot[ , 6 ] %in% c( "-1", -1, "-" ), ] <- rev( cov_df[ annot[ , 6 ] %in% c( "-1", -1, "-" ), ] )
  }

  #~ colnames( cov_df ) <- c( -999:1000 )
  start_pos=( 1 - off_set )
  colnames( cov_df ) <- seq( from=start_pos, to=( start_pos + ncol( cov_df ) - 1 )  ) #~  c( ( 1 - ceiling( annot[ 1, 3 ] - annot[ 1, 2 ] ) / 2 ) : ( floor( annot[ 1, 3 ] - annot[ 1, 2 ] ) / 2 ) )

  # harmonize the scale among the exons
  if ( harmonization ) {
    idx <- which( rowSums( cov_df ) != 0 )
    cov_df[ idx, ] <- as.data.frame( t( apply( cov_df[ idx, ], 1, function( xxx ) { xxx / max( xxx ) } ) ) )
  }

  return( cov_df )
}
