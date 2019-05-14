#!/usr/bin/Rscript

require( 'ggplot2' )

#### build measure metagene
metagene_measure_builder <- function( filename, melt_df, ylims=NULL, diff=FALSE, colour=NULL, linetype=NULL, color_pallette=NULL ) {
  if ( is.null( colour ) ) {
    colour <- melt_df$replicate
  }
  leg_nrow <- 2 #( length( levels( factor( colour ) ) ) + 1 ) %/% 2


  if ( is.null( linetype ) ) {
    linetype <- melt_df$condition
  }

  fig_aes <- aes( x = X1, y = value, colour = colour, linetype = linetype )
  if ( length( levels( as.factor( linetype ) ) ) == 1 ) {
    fig_aes$linetype <- NULL
  }
  figure <- ggplot( data = melt_df, fig_aes )

  # if ( diff ) {
  #   figure <- ggplot( data = melt_df, aes( x = X1, y = value, colour = colour ) )
  # } else {
  #   figure <- ggplot( data = melt_df, aes( x = X1, y = value, colour = colour, linetype = linetype ) )
  # }

  # png( filename = filename, width = 1920, height = 1200, pointsize = 48 )
  extrem_values <- c(
    min( melt_df[ , "value" ] ),
    max( melt_df[ , "value" ] )
  )

  if ( length( ylims ) != 2 ) {
    ylims <- extrem_values
  } else {
    if ( is.na( ylims[ 1 ] ) ) {
      ylims[ 1 ] <- extrem_values[ 1 ]
    }
    if ( is.na( ylims[ 2 ] ) ) {
      ylims[ 2 ] <- extrem_values[ 2 ]
    }
  }


  message( '>>> metaplot' )
  message( paste0( 'min value in plotted data: ', extrem_values[ 1 ] ) )
  message( paste0( 'max value in plotted data: ', extrem_values[ 2 ] ) )


  # ylims=c( 0, 100 )
  figure <- figure + theme(text = element_text(size=48)) + coord_cartesian( ylim=ylims ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))
  # + scale_y_continuous( trans='log10' )

  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_colour_manual( values=color_pallette )
  }


  message( filename )
  png( filename = filename )
  print(
    figure + geom_line(size=1.5)
  )

  dev.off()
}
