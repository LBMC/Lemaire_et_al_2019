#!/usr/bin/Rscript

cov_extr_fc_dir <- paste( dirname(sys.frame(1)$ofile), 'metagene_coverage_functions', sep='/' )
source( paste( cov_extr_fc_dir, 'stranded_coverages.r', sep='/' ) )
source( paste( cov_extr_fc_dir, 'stranded_diffCov.r', sep='/' ) )


## extract the coverage from a batch of samples for one list of annotations
cov_extr_fc <- function( bw_dir, comp_pair, seqDepth_df, harmonization=FALSE ) {
  cov_df_list <- list()
  for ( xxx in rep_levels ) {
    if ( length( comp_pair ) == 1 ) {
      bw_pattern <- paste( comp_pair, 'n', xxx, ".*bw$", sep='' )
      bw_name <- list.files( path = bw_dir, pattern = bw_pattern )

      bw_file <- paste( bw_dir, bw_name, sep='/' )

      cov_df_list[[ xxx ]] <- stranded_coverages( bw_file, annot, harmonization=harmonization )
    } else if ( length( comp_pair ) == 2 ) {
      bw_file_list <- list()
      normFactor_list <- list()
      for ( cond in comp_pair ) {
        bw_pattern <- paste( cond, 'n', xxx, ".*bw$", sep='' )
        bw_name <- list.files( path = bw_dir, pattern = bw_pattern )

        bw_file_list[[ cond ]] <- paste( bw_dir, bw_name, sep='/' )
      }

      cov_df_list[[ xxx ]] <- stranded_diffCov( bw_file_list=bw_file_list, annot=annot, comp_pair=comp_pair, harmonization=TRUE )
    }
  }
  return( cov_df_list )
}
