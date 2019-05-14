#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

exon_length_repart_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( exon_length_repart_plots_dir, 'fig_functions_source.r', sep='/' ) )
source( paste( exon_length_repart_plots_dir, 'fig_format.r', sep='/' ) )
fig_functions_source( exon_length_repart_plots_dir )

exon_length_repart_plots <- function ( final_data_list, nom_vec=c( 'query' ), nom_query=c( 'query' ) ,out_dir, sign, color_pallette=NULL ) {
  dir.create( out_dir, showWarnings=FALSE, recursive=TRUE )

  col <- c( 'EXON_LENGTH_BEFORE', 'EXON_LENGTH', 'EXON_LENGTH_AFTER' )
  if ( all( col %in% colnames( final_data_list[[ nom_query[ 1 ] ]] ) ) ) {
    for ( xxx in names( final_data_list ) ) {
      for ( col_name in c( 'EXON_LENGTH_AFTER', 'EXON_LENGTH_BEFORE' ) ) {
        if ( ! col_name %in% colnames( final_data_list[[ xxx ]] ) ) {
          final_data_list[[ xxx ]][ col_name ] <- final_data_list[[ xxx ]][ 'EXON_LENGTH' ]
        }
      }
    }

    for ( col_name in col ) {
      data_tab <- ggplot_repart_tab( final_data_list, col=col_name, nom_vec=nom_vec )
      leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2

      fig <- ggplot( data=data_tab, aes( x=value, y=repart, color=liste ) ) + scale_x_continuous( trans='log10' ) + ylab( 'fraction') + xlab( 'value' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
      #~ coord_cartesian( xlim=c( 0, 10000 ) )

      # set colors if defined
      if ( ! is.null( color_pallette ) ) {
        fig <- fig + scale_colour_manual( values=color_pallette )
      }

      png( fig_name( col_name, out_dir, sign, 'repart'), width=fig_width( col ) )
      print(
        fig + geom_point() + geom_line()
      )
      bouh <- dev.off()
    }

  }

}
