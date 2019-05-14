#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

FarLine_psi_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( FarLine_psi_plots_dir, 'fig_functions_source.r', sep='/' ) )
fig_functions_source( FarLine_psi_plots_dir )

FarLine_psi_plots <- function ( list_psi_tab, nom_vec=c( 'query' ), out_dir, sign='query', suffix='PSI', draw_quantiles=c( 0.25, 0.5, 0.75 ), color_pallette=NULL ) {

  data_tab <- ggplot_data_tab( list_psi_tab, col=NULL, nom_vec )
  col <- data_tab$variable[ !duplicated( data_tab$variable ) ]
  leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2


  # build the 2D + density plots for each list
  sub_fig_name_root <- unlist( strsplit( fig_name( suffix, out_dir=out_dir, nom=sign, '' ), split='.png', fixed=TRUE ) )[ 1 ]
  for ( annot_name in nom_vec ) {
    sub_fig_name <- paste0( sub_fig_name_root, annot_name, '.png' )
    axes_names <- colnames( list_psi_tab[[ annot_name ]] )
    sub_fig <- ggplot( list_psi_tab[[ annot_name ]], aes_string( x=axes_names[ 2 ], y=axes_names[ 1 ] ) )
    sub_fig <- sub_fig + xlim( 0, 1 ) + ylim( 0, 1 ) + ylab( paste0( axes_names[ 1 ], ' PSI median' ) ) + xlab( paste0( axes_names[ 2 ], ' PSI median' ) ) + theme_bw( base_size = 16 ) + theme(legend.position="none")

    png( sub_fig_name )
    print(
      sub_fig + geom_abline( intercept=0, slope=1 ) + geom_point() + geom_density2d()
    )
    bouh <- dev.off()
  }


  ## build the boxplots/violin plots
  fig <- ggplot( data=data_tab, aes( x=variable, y=value, fill=liste ) ) + ylim( 0, 1 ) + ylab( 'Percentage Splicing Inclusion (PSI)') + xlab( 'condition' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = leg_nrow))

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    fig <- fig + scale_fill_manual( values=color_pallette )
  }

png( fig_name( suffix, out_dir=out_dir, nom=sign, 'violin'), width=fig_width( col ) )
  print(
    fig + geom_violin( draw_quantiles=draw_quantiles )
    )
  bouh <- dev.off()

  png( fig_name( suffix, out_dir=out_dir, nom=sign, 'boxplot'), width=fig_width( col ) )
  print(
    fig + geom_boxplot()
    )
  bouh <- dev.off()
}
