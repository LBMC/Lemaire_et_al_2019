#!/usr/bin/Rscript

require( 'ggplot2' )

#### build summary metagene
metagene_summary_builder <- function( filename, melt_df, ylims=NULL ) {
  figure <- ggplot( data = melt_df, aes( x = X1, y = value, colour = measure, linetype = replicate ) )
  png( filename = filename, width = 1920, height = 1200, pointsize = 48 )

  if ( length( ylims ) != 2 ) {
    ylims <- c( min( melt_df[ , "value" ] ), max( melt_df[ , "value" ] ) )
  }

  print(
    figure + geom_line(size=1.5) + theme(text = element_text(size=48)) + ylim( ylims )
  )

  dev.off()
}
