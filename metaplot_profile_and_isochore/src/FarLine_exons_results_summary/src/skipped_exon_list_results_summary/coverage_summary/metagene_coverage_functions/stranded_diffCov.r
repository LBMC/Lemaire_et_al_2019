#!/usr/bin/Rscript

require( 'bigWig' )

stranded_diffCov_dir <- paste( dirname(sys.frame(1)$ofile), 'metagene_coverage_functions', sep='/' )
source( paste( stranded_diffCov_dir, 'stranded_coverages.r', sep='/' ) )

## differential coverage between two conditions
stranded_diffCov <- function( bw_file_list, annot, off_set=1, normFactor_list=list(), comp_pair, harmonization = FALSE ) {
  cov_df_list <- list()
  for ( cond in comp_pair ) {
    if ( is.null( normFactor_list[[ cond ]] ) ) {
      normFactor_list[[ cond ]] <- 1
    }

    cov_df_list[[ cond ]] <- stranded_coverages( bw_file_list[[ cond ]], annot, off_set=off_set, normFactor=normFactor_list[[ cond ]] )
  }

  # diff_cov_df <- data.frame( cov_df_list[[ comp_pair[ 1 ] ]] - cov_df_list[[ comp_pair[ 2 ] ]], check.names=FALSE )
  # message( c( "--- operation: ", comp_pair[ 2 ], '-' , comp_pair[ 1 ] ))
  diff_cov_df <- cov_df_list[[ comp_pair[ 2 ] ]] - cov_df_list[[ comp_pair[ 1 ] ]]
  if ( dim( annot )[ 2 ] >= 6 ) {
    diff_cov_df[ annot[ , 6 ] %in% c( "-1", -1, "-" ), ] <- rev( diff_cov_df[ annot[ , 6 ] %in% c( "-1", -1, "-" ), ] )
  }

  # harmonize the scale among the exons
  if ( harmonization ) {
    idx <- which( rowSums( abs( diff_cov_df ) ) != 0 )
    diff_cov_df[ idx, ] <- as.data.frame( t( apply( diff_cov_df[ idx, ], 1, function( xxx ) { xxx / max( abs( xxx ) ) } ) ) )
  }

  return( diff_cov_df )
}
