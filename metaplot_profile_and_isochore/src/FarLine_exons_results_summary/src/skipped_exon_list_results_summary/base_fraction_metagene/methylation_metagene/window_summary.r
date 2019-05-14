#!/usr/bin/Rscript

## average the methylation rate in a scanning window
window_summary <- function( mat, window=1 ) {
  if ( window > 1 ) {
    for ( xxx in 1:( ncol( mat ) - ( window - 1 ) ) ) {
      mat[ , xxx ] <- rowMeans( mat[ , xxx:( xxx + window - 1 ) ], na.rm=TRUE )
    }

    shift_start <- floor( window / 2.0 )
    shift_end <- ( ceiling( window / 2.0 ) - 1 )
    colonnes <- colnames( mat )[ ( 1 + shift_start ):( ncol( mat ) - shift_end ) ]

    mat <- mat[ , 1:( ncol( mat ) - ( window - 1 ) ) ]
    colnames( mat ) <- colonnes
  }

  return( mat )

}
