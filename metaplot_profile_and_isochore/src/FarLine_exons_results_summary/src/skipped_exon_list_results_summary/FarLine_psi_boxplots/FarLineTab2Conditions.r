#!/usr/bin/Rscript

FarLineTab2Conditions <- function ( FarLine_all_data_tab ) {
  psi_columns <- colnames( FarLine_all_data_tab )[ grep( pattern='_psi$', x=colnames( FarLine_all_data_tab ) ) ]
  conditions <- as.data.frame( do.call( rbind, strsplit( psi_columns, split='_' ) ), stringsAsFactors=FALSE )
  if ( ncol( conditions ) == 3 ) {
    conditions <- conditions[ , 1 ]
  } else {
    conditions <- apply( conditions[ , 1:( ncol( conditions ) - 2 ) ], 1, paste, collapse='_' )
  }
  conditions <- conditions[ !duplicated( conditions ) ]
  return( conditions )
}
