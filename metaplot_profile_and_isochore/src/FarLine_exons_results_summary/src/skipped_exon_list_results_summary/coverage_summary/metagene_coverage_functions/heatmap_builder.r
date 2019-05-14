#!/usr/bin/Rscript

#### build heatmap
heatmap_builder <- function( filename, cov_df, scale="none", Rowv=NA, Colv=NA, cexCol=0.2, col=heat.colors( 256 ), ymax=NULL, sum.order=FALSE, ... ) {
    # message( filename )

    if ( sum.order ) {
      sum_order=order( apply( cov_df, 1, sum, na.rm=TRUE ) )
      cov_df <- cov_df[ sum_order, ]
    }

    if ( ! is.null( ymax ) ) {
      cov_df[ cov_df > ymax ] <- ymax
    }

    extrem_values <- c(
      min( as.matrix( cov_df ) ),
      max( as.matrix( cov_df ) )
    )

    message( '>>> heatmap')
    message( paste0( 'min value in plotted data: ', extrem_values[ 1 ] ) )
    message( paste0( 'max value in plotted data: ', extrem_values[ 2 ] ) )

    message( filename )
    png( filename = filename, width = 1200, height = 1200, pointsize = 48 )
    heatmap( as.matrix( cov_df ), cexCol=cexCol, Colv=Colv, Rowv=Rowv, scale=scale, col=col, ... )
    dev.off()
}
