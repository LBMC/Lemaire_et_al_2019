#!/usr/bin/Rscript

require( ggplot2 )

base_fraction_repart_plot_dir <- dirname(sys.frame(1)$ofile)
source( paste( base_fraction_repart_plot_dir, 'base_fraction_metagene_fig_name.r', sep='/' ) )
source( paste( base_fraction_repart_plot_dir, 'base_colors.r', sep='/' ) )

## build the curves of metagenes
base_fraction_repart_plot_fc <- function( list_fraction_mat, vec_annot_name, color_pallette=NULL, base_vec, skew_measure=NULL, prefix='', entropy=FALSE, struct_info=FALSE, ent_base_vec=NULL, regexpr_fig_name=NULL, presence.only=FALSE ) {
  list_ratio_vec <- list()

  # compute the name of the figure (and create the output folder)
  metagene_fig_name <- base_fraction_metagene_fig_name_fc( base_vec, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )
  metagene_dir <- dirname( metagene_fig_name )
  repart_dir <- paste0( metagene_dir, '/reparts' )
  dir.create( repart_dir, showWarnings = FALSE, recursive = TRUE )

  ## repartition for each annotation list
  for ( annot_name in vec_annot_name ) {
    # out_prefix <- paste0( prefix, '_', annot_name )
    out_prefix <- paste0( repart_dir, '/', unlist( strsplit( basename( metagene_fig_name ), split='.', fixed=TRUE ) )[ 1 ], '_', annot_name )
    fig_name <- paste0( out_prefix, '.png' )
    # message( fig_name )

    # compute the sum of signal
    sum_vec <- apply( list_fraction_mat[[ annot_name ]], 1, sum, na.rm=TRUE )
    cum_fract <- order( order( sum_vec ) ) / length( sum_vec )
    ratio_vec <- sum_vec / ncol( list_fraction_mat[[ annot_name ]] )
    list_ratio_vec[[ annot_name ]] <- data.frame( ratio_vec=ratio_vec, cum_fract=cum_fract )

    # build the repart
    png( fig_name )
    plot( ratio_vec, cum_fract, pch=16 )
    bouh <- dev.off()
  }

  ## plot with all the repartitions
  # reshape the list of metagenes in a data.frame
  df_ratio <- melt( list_ratio_vec, id.vars=c( 'cum_fract', 'ratio_vec' ) ) #listMeta2DF( list_ratio_vec )
  df_ratio$L1 <- factor( df_ratio$L1, levels=vec_annot_name )
  leg_nrow <- ( length( levels( factor( df_ratio$L1 ) ) ) + 1 ) %/% 2

  # build the plot
  figure <- ggplot( data=df_ratio, aes( x=ratio_vec, y=cum_fract, color=L1 ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_colour_manual( values=color_pallette )
  }

  # set the name of the figure
  fig_name <- paste0( repart_dir, '/', basename( metagene_fig_name ) )
  # message( fig_name )


  # write the figure
  dir.create( dirname( fig_name ),showWarnings = FALSE, recursive = TRUE )
  png( fig_name )
  print(
    figure + geom_line( size=1.5 ) + geom_point( size=1.5 ) + xlab( 'signal value' ) + ylab( 'cummulative fraction' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom")  + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))
    # + ylim( 0, 1 )
  )
  bouh <- dev.off()

}
