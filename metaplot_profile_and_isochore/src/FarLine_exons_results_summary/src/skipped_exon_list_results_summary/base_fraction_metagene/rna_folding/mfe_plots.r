#!/usr/bin/Rscript

mfe_plots <- function( mfe_df_list, names_vec, fig_prefix, color_pallette=NULL, draw_quantiles=c( 0.25, 0.5, 0.75 ) ) {
  gg_data <- melt( mfe_df_list, id.vars=c( 'id', 'MFE', 'fraction' ) )
  gg_data$liste <- factor( gg_data$liste, levels=names_vec )
  leg_nrow <- ( length( levels( factor( gg_data$liste ) ) ) + 1 ) %/% 2

  fig <- ggplot( gg_data )
  fig <- fig + theme_bw( base_size = 16 ) + theme(legend.position="bottom")  + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))

  ## repartitions
  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    fig <- fig + scale_colour_manual( values=color_pallette )
  }

  figure <- fig + geom_point( mapping=aes( x=MFE, y=fraction, col=liste ) ) + geom_line( mapping=aes( x=MFE, y=fraction, col=liste ) ) + xlab( 'MFE' ) + ylab( 'fraction' )

  fig_name <- paste0( fig_prefix, '_repart.png')
  # message( fig_name )

  png( fig_name, width=800, height=800 )
  print(
    figure
  )
  bouh <- dev.off()
  ####


  ## statistical tests


  ## boxplots and violin plots
  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    fig <- fig + scale_fill_manual( values=color_pallette )
  }

  boxplot_fig <- fig + geom_boxplot( mapping=aes( x=liste, y=MFE, fill=liste ) ) + xlab( 'exon liste' ) + ylab( 'MFE' )

  fig_name <- paste0( fig_prefix, '_boxplot.png')
  # message( fig_name )

  png( fig_name, width=800, height=800 )
  print(
    boxplot_fig
  )
  bouh <- dev.off()


  violin_fig <- fig + geom_violin( mapping=aes( x=liste, y=MFE, fill=liste ), draw_quantiles=draw_quantiles ) + xlab( 'exon liste' ) + ylab( 'MFE' )

  fig_name <- paste0( fig_prefix, '_violin.png')
  # message( fig_name )

  png( fig_name, width=800, height=800 )
  print(
    violin_fig
  )
  bouh <- dev.off()



  ####
}
