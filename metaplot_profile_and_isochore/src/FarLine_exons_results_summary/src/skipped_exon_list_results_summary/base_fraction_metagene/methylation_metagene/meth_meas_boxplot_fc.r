#!/usr/bin/Rscipt

require( ggplot2 )
require( reshape2 )

meth_meas_boxplot_fc <- function( all_list_meth_meas, out_prefix, vec_annot_name, color_pallette=NULL ) {

  ggdata <- melt( all_list_meth_meas )
  ggdata[ 'group' ] <- paste( ggdata$L1, ggdata$L3, ggdata$L2, sep='_' )
  ggdata$L3 <- factor( ggdata$L3, levels=vec_annot_name )

  colnames( ggdata ) <- c( 'value', 'index', 'liste', 'replicate', 'condition', 'group' )

  ylim_basal <- 0
  ylabel <- 'methylation rate'

  for ( cond in unique( ggdata$condition ) ) {
    # build the figure
    f1 <- ggplot( ggdata[ ggdata$condition == cond, ], aes( x=group, y=value, group=group, col=liste ) ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + labs( y=ylabel ) + ylim( ylim_basal, 1 )

    # set colors if defined
    if ( ! is.null( color_pallette ) ) {
      f1 <- f1 + scale_fill_manual( values=color_pallette ) + scale_colour_manual( values=color_pallette )
    }

    # add plots to the figures
    f1 <- f1 + geom_violin( adjust=1/2, scale='width' ) + stat_summary(fun.data="median_hilow", geom="pointrange", color="black")

    # write the figure
    fig_name <- paste0( out_prefix, '_', cond, '_mCpG_dist.png' )
    png( fig_name, width=600, height=600 )
    print(
      f1
    )
    bouh <- dev.off()
  }
}
