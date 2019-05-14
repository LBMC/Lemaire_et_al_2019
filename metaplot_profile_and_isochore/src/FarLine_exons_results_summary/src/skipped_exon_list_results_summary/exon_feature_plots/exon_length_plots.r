#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

exon_length_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( exon_length_plots_dir, 'fig_functions_source.r', sep='/' ) )
source( paste( exon_length_plots_dir, 'fig_format.r', sep='/' ) )
fig_functions_source( exon_length_plots_dir )

exon_length_plots <- function ( final_data_list, nom_vec=c( 'query' ), nom_query=c( 'query' ), out_dir, sign, draw_quantiles=c( 0.25, 0.5, 0.75 ), color_pallette=NULL ) {
  col <- c( 'EXON_LENGTH_BEFORE', 'EXON_LENGTH', 'EXON_LENGTH_AFTER' )
  if ( all( col %in% colnames( final_data_list[[ nom_query[ 1 ] ]] ) ) ) {
    for ( xxx in names( final_data_list ) ) {
      for ( col_name in c( 'EXON_LENGTH_AFTER', 'EXON_LENGTH_BEFORE' ) ) {
        if ( ! col_name %in% colnames( final_data_list[[ xxx ]] ) ) {
          final_data_list[[ xxx ]][ col_name ] <- final_data_list[[ xxx ]][ 'EXON_LENGTH' ]
        }
      }
    }

    data_tab <- ggplot_data_tab( final_data_list, col, nom_vec )
    leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2

    fig <- ggplot( data=data_tab, aes( x=variable, y=value, fill=liste ) ) + scale_y_continuous( trans='log10' ) + ylab( 'value') + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
    #~_coord_trans( y='log10' )
    #~ coord_cartesian( ylim=c( 0, 400 ) )

    # set colors if defined
    if ( ! is.null( color_pallette ) ) {
      fig <- fig + scale_fill_manual( values=color_pallette )
    }

    png( fig_name( 'exon_length', out_dir, sign, 'violin'), width=fig_width( col ) )
    print(
      fig + geom_violin( draw_quantiles=draw_quantiles )
      )
    dev.off()

    png( fig_name( 'exon_length', out_dir, sign, 'boxplot'), width=fig_width( col ) )
    print(
      fig + geom_boxplot()
      )
    dev.off()
  }
}
