#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

diffDonor_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( diffDonor_plots_dir, 'fig_functions_source.r', sep='/' ) )
source( paste( diffDonor_plots_dir, 'fig_format.r', sep='/' ) )
fig_functions_source( diffDonor_plots_dir )

diffDonor_plots <- function ( final_data_list, nom_vec=c( 'query' ), out_dir, sign, draw_quantiles=c( 0.25, 0.5, 0.75 ), color_pallette=NULL ) {
  col <- c( 'STRENGTH_DONOR_XvsBEFORE', 'STRENGTH_DONOR_XvsAFTER' )
  if ( all( col %in% colnames( final_data_list[[ 1 ]] ) ) ) {

    data_tab <- ggplot_data_tab( final_data_list, col, nom_vec )
    leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2

    fig <- ggplot( data=data_tab, aes( x=variable, y=value, fill=liste ) ) + scale_y_continuous( trans='identity', limits=c( 0, NA ) ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))

    # set colors if defined
    if ( ! is.null( color_pallette ) ) {
      fig <- fig + scale_fill_manual( values=color_pallette )
    }

    png( fig_name( 'force_diffDonor', out_dir, sign, 'violin'), width=fig_width( col ) )
    print(
      fig + geom_violin( draw_quantiles=draw_quantiles )
      )
    dev.off()

    png( fig_name( 'force_diffDonor', out_dir, sign, 'boxplot'), width=fig_width( col ) )
    print(
      fig + geom_boxplot()
      )
    dev.off()
  }
}
