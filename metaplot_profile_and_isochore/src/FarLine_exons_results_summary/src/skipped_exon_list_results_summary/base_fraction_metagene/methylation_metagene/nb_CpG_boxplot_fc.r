#!/usr/bin/Rscipt

require( ggplot2 )
require( reshape2 )

nb_CpG_boxplot_fc <- function( list_nb_CpG, out_prefix, vec_annot_name, color_pallette=NULL ) {

  ggdata <- melt( list_nb_CpG )
  ggdata$L1 <- factor( ggdata$L1, levels=vec_annot_name )

  # build the figure
  f1 <- ggplot( ggdata, aes( x=L1, y=value, col=L1 ) ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + labs( x='liste', y='fraction CpG' )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    f1 <- f1 + scale_fill_manual( values=color_pallette ) + scale_colour_manual( values=color_pallette )
  }

  # add plots to the figures
  f1 <- f1 + geom_violin( adjust=1/2, scale='width' ) + stat_summary(fun.data="median_hilow", geom="pointrange", color="black")

  # write the figure
  fig_name <- paste0( out_prefix, '_nbCpG.png' )
  png( fig_name, width=600, height=600 )
  print(
    f1
  )
  bouh <- dev.off()
}
