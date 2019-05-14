#!/usr/bin/Rscript

require( 'ggplot2', quietly = TRUE )
require( 'reshape2', quietly = TRUE )


base_fraction_metagene_plot_dir <- dirname(sys.frame(1)$ofile)
source( paste( base_fraction_metagene_plot_dir, 'base_fraction_metagene_fig_name.r', sep='/' ) )

## build the curves of metagenes
base_fraction_metagene_plot_fc <- function( list_metagene, vec_annot_name, base_vec, color_pallette=NULL, skew_measure=NULL, prefix='./', entropy=FALSE, struct_info=FALSE, ent_base_vec=NULL, regexpr_fig_name=NULL, presence.only=FALSE, ylim_arg=NULL ) {
  # reshape the list of metagenes in a data.frame
  df_metagene <-listMeta2DF( list_metagene )
  df_metagene$L1 <- factor( df_metagene$L1, levels=vec_annot_name )
  leg_nrow <- ( length( levels( factor( df_metagene$L1 ) ) ) + 1 ) %/% 2

  # build the plot
  figure <- ggplot( data=df_metagene, aes( x=position, y=value, color=L1 ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_colour_manual( values=color_pallette )
  }

  # set Y-axis range if defined
  if ( ! is.null( ylim_arg ) ) {
    figure <- figure + ylim( ylim_arg[1], ylim_arg[2] )
  }

  # set the name of the figure
  # message( fig_name )
  fig_name <- base_fraction_metagene_fig_name_fc( base_vec, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )

  # write the figure
  dir.create( dirname( fig_name ),showWarnings = FALSE, recursive = TRUE )
  png( fig_name )
  print(
    figure + geom_line( size=1.5 ) + xlab( 'relative position to splicing site (bp)' ) + ylab( 'fraction' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom")  + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))
    # + ylim( 0, 1 )
  )
  bouh <- dev.off()
}
