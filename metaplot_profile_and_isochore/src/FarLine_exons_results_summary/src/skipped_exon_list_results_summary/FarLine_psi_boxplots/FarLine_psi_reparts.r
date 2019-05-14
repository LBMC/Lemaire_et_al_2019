#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

FarLine_psi_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( FarLine_psi_plots_dir, 'fig_functions_source.r', sep='/' ) )
fig_functions_source( FarLine_psi_plots_dir )

FarLine_psi_reparts <- function ( list_psi_repart_tab, nom_vec=c( 'query' ), out_dir, sign='query', color_pallette=NULL ) {
  data_tab <- melt( list_psi_repart_tab, id.vars=c( 'psi', 'repart' ) ) #ggplot_data_tab( list_psi_repart_tab, col=NULL, nom_vec )
  colnames( data_tab ) <- c( 'psi', 'repart', 'cond', 'liste' )
  data_tab$liste <- factor( data_tab$liste, levels=nom_vec )
  col <- data_tab$variable[ !duplicated( data_tab$variable ) ]
  leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2


  for ( cond in unique( data_tab$cond ) ) {
    ## build the boxplots/violin plots
    fig <- ggplot( data=data_tab[ data_tab$cond == cond, ], aes( x=psi, y=repart, color=liste ) ) + ylim( 0, 1 ) + xlab( 'max( PSI )') + ylab( 'fraction' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = leg_nrow))

    # set colors if defined
    if ( ! is.null( color_pallette ) ) {
      fig <- fig + scale_color_manual( values=color_pallette )
    }

    png( fig_name( cond, out_dir=out_dir, nom=sign, 'PSIreparts') ) #, width=fig_width( col ) )
    print(
      fig + geom_point() + geom_line()
    )
    bouh <- dev.off()
  }
}
