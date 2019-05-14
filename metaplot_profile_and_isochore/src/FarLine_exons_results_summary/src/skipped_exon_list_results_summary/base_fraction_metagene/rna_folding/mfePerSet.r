#!/usr/bin/Rscript

mfePerSet_dir <- dirname(sys.frame(1)$ofile)
source( paste( mfePerSet_dir, 'rnafold2MFEtable.r', sep='/' ) )

mfePerSet <- function( RNAfold_res_vec, names_vec=NULL ) {
  mfe_df_list <- list(); attr( mfe_df_list, "varname") <- "liste"

  if ( is.null( names_vec ) ) {
    names_vec <- as.character( 1:length( RNAfold_res_vec ) )
  }

  for ( filepath_idx in 1:length( RNAfold_res_vec ) ) {
    filepath <- RNAfold_res_vec[ filepath_idx ]
    name_set <- names_vec[ filepath_idx ]

    mfe_df_list[[ name_set ]] <- rnafold2MFEtable( filepath )

    ## compute the repartition
    mfe_df_list[[ name_set ]][ 'fraction' ] <- order( order( mfe_df_list[[ name_set ]][ , 'MFE' ] ) ) / nrow( mfe_df_list[[ name_set ]] )

  }

  return( mfe_df_list )
}
