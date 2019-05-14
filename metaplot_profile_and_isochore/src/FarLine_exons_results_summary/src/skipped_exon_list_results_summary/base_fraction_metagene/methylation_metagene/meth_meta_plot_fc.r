#!/usr/bin/Rscipt

require( ggplot2 )
require( reshape2 )

meth_meta_plot_fc <- function( list_meth_meta_list, out_prefix, vec_annot_name, color_pallette=NULL, diff.mode=FALSE ) {
  # reshape the data
  # list_meth_meta_list <- lapply( list_meth_meta_list, function( xxx ) {
  #   lapply( xxx, function( yyy ) {
  #     data.frame( rel_pos=as.numeric( names( yyy ) ), value=yyy )
  #   })
  # })

  ggdata <- melt( list_meth_meta_list, id.vars=c( 'rel_pos', 'value', 'CpG_counts' ) )
  ggdata$L1 <- factor( ggdata$L1, levels=vec_annot_name )

  ylim_basal <- 0
  ylabel <- 'methylation rate'
  if ( diff.mode ) {
    ylim_basal <- -1
    ylaber <- 'diff. methylation rate'
  }

  # build the figure
  f1 <- ggplot( ggdata, aes( x=rel_pos, y=value, color=L1 ) ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + labs( y=ylabel ) + ylim( ylim_basal, 1 )
  f2 <- ggplot( ggdata, aes( x=rel_pos, y=CpG_counts, color=L1 ) ) + theme_bw( base_size = 16 ) + theme(legend.position="none") + labs( x='relative position\nto splicing site', y='fraction exon with CpG' )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    f1 <- f1 + scale_colour_manual( values=color_pallette )
    f2 <- f2 + scale_colour_manual( values=color_pallette )
  }

  # add plots to the figures
  f1 <- f1 + geom_line( size=1.5 )
  f2 <- f2 + geom_line( size=1.5 )

  # write the figure
  fig_name <- paste0( out_prefix, '_mCpG.png' )
  png( fig_name, width=800, height=800 )
  grid.arrange( f1, f2, heights=c( 4/5, 1/5 ) )
  bouh <- dev.off()
}
