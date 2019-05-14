#!/usr/bin/Rscript

require( ggplot2 )
require( purrr )
coverage_repart_plot_dir <- dirname(sys.frame(1)$ofile)

## build the curves of metagenes
coverage_repart_plot_fc <- function( list_vec, vec_annot_name, color_pallette=NULL, filename='./test.png' ) {
  list_ratio_vec <- list()

  fig_dir <- dirname( filename )
  dir.create( fig_dir, showWarnings = FALSE, recursive = TRUE )

  ## plot with all the repartitions
  # reshape the list of metagenes in a data.frame
  list_vec <- modify_depth( list_vec, 3, function( xxx ) {
    data.frame( max_val=xxx, fraction=order( order( xxx ) ) / length( xxx ) )
  } )
  df_data <- melt( list_vec, id.vars=c( 'max_val', 'fraction' ) ) #listMeta2DF( list_ratio_vec )
  df_data$L1 <- factor( df_data$L1, levels=vec_annot_name )
  df_data[ 'combi' ] <- paste( df_data$L1, df_data$L2, df_data$L3, sep='_' )
  leg_nrow <- ( length( levels( factor( df_data$L1 ) ) ) + 1 ) %/% 2

  # build the plot
  figure <- ggplot( data=df_data, aes( x=max_val, y=fraction, color=L1, group=df_data$combi ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_colour_manual( values=color_pallette )
  }


  # write the figure
  png( filename )
  message( filename )
  print(
    figure + geom_line( size=1.5 ) + geom_point( size=1.5 ) + xlab( 'max value' ) + scale_x_continuous( trans='log10' ) + ylab( 'cummulative fraction' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom")  + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))
    # + ylim( 0, 1 )
  )
  bouh <- dev.off()

}
