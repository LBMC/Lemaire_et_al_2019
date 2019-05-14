#!/usr/bin/Rscript

FarLine_psi_summary_dir <- dirname(sys.frame(1)$ofile)

## Some function
source( paste( FarLine_psi_summary_dir, 'FarLineTab2Conditions.r', sep='/' ) )


FarLine_psi_summary <- function ( FarLine_all_data_tab, reduc_fc=median, ... ) {
  conditions <- FarLineTab2Conditions( FarLine_all_data_tab )

  out_tab <- as.data.frame( matrix( NA, nrow=nrow( FarLine_all_data_tab ), ncol=length( conditions ), dimnames=list( NULL, conditions ) ) )
  for ( cond in conditions ) {
    pattern <- paste( '^', cond, '.*B[0-9]([0-9]*)?s_psi$', sep='' )
    query_col <- colnames( FarLine_all_data_tab )[ grep( pattern=pattern, x=colnames( FarLine_all_data_tab ) ) ]
    out_tab[ cond ] <- apply( FarLine_all_data_tab[ , query_col ], 1, reduc_fc, ... )
  }
  return( out_tab )
}
