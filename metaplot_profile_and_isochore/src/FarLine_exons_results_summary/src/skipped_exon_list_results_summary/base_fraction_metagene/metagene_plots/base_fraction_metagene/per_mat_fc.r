#!/usr/bin/Rscript

###########
per_mat_fc <- function( seq_mat ) {
  per_mat <- list()
  # per_mat <- matrix( NA, nrow=nrow( seq_mat ), ncol=ncol( seq_mat ) )

  for ( annot in 1:nrow( seq_mat ) ) {
    per_mat[[ annot ]] <- acf( seq_mat[ annot, ], plot=FALSE, na.action=na.pass )[[ 'acf' ]][ , 1, 1 ]
  }
  per_mat <- do.call( rbind, per_mat )

  return( list( 'per_mat'=per_mat, 'acf_struct'=acf( seq_mat[ 1, ], plot=FALSE, na.action=na.pass ) ) )
}
