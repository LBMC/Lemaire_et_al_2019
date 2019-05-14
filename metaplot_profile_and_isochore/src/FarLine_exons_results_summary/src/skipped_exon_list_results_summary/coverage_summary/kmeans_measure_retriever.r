#!/usr/bin/Rscript

kmeans_measure_retriever <- function ( df_expdesign, seqDepth_df, annot, repl, comp_pair, cov_harm=FALSE, diff_harm=TRUE ) {
  if ( length( comp_pair ) == 1 ) {
    cond <- comp_pair
    bw_file <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]

    # retrieve the matrix of measure
    cov_df <- stranded_coverages( bw_file, annot, seqDepth_df[ seqDepth_df$condition==comp_pair & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=cov_harm )

  } else if ( length( comp_pair ) == 2 ) {
    bw_file_list <- list()
    normFactor_list <- list()

    for ( cond in comp_pair ) {
      # message( cond )
      bw_file_list[[ cond ]] <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]
      normFactor_list[[ cond ]] <- seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ]
    }

    # retrieve the matrix of measure
    cov_df <- stranded_diffCov( bw_file_list, annot, normFactor_list, comp_pair, harmonization=diff_harm )
  }

  return( cov_df )
}
