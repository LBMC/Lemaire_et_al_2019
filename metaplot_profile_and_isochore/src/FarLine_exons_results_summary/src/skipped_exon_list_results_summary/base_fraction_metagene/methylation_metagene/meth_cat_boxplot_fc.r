#!/usr/bin/Rscipt

require( ggplot2 )
require( reshape2 )

meth_cat_boxplot_fc <- function( all_list_meth_cat, out_prefix, vec_annot_name, color_pallette=NULL ) {

  ggdata <- melt( all_list_meth_cat )
  colnames( ggdata ) <- c( 'value', 'category', 'index', 'liste', 'replicate', 'condition' )
  ggdata[ 'group' ] <- paste( ggdata$condition, ggdata$liste, ggdata$replicate, sep='_' )
  ggdata$liste <- factor( ggdata$liste, levels=vec_annot_name )
  ggdata$category <- factor( ggdata$category, levels=c( 'Meth', 'Inter', 'Raw' ) )


  for ( cond in unique( ggdata$condition ) ) {
    # build the figure
    f1 <- ggplot( ggdata[ ggdata$condition == cond, ], aes( x=group, y=value, group=group, col=liste ) ) + facet_wrap( ~ category, ncol=1 ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank())

    # set colors if defined
    if ( ! is.null( color_pallette ) ) {
      f1 <- f1 + scale_fill_manual( values=color_pallette ) + scale_colour_manual( values=color_pallette )
    }

    # add plots to the figures
    f1 <- f1 + geom_violin( adjust=1/2, scale='width' ) + stat_summary(fun.data="median_hilow", geom="pointrange", color="black") #, draw_quantiles=c( 0.25, 0.5, 0.75 )

    # write the figure
    fig_name <- paste0( out_prefix, '_', cond, '_CpG-cat.png' )
    png( fig_name, width=600, height=600 )
    print(
      f1
    )
    bouh <- dev.off()
  }
}
