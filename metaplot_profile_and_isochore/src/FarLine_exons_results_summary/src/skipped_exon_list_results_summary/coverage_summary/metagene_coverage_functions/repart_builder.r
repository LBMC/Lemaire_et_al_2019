#!/usr/bin/Rscript

#### build heatmap
repart_builder <- function( filename, val_vec, ... ) {
    # message( filename )

    message( '>>> repartMax')
    message( 'summary: ' )
    sum_res <- summary( val_vec )
    message( paste( names( sum_res ), collapse='\t' ) )
    message( paste( ( sum_res ), collapse='\t' ) )

    png( filename )
    plot( max_vec[ order( val_vec ) ], ( 1:length( val_vec ) ) / length( val_vec ), pch=16, ... )
    bouh <- dev.off()

}

