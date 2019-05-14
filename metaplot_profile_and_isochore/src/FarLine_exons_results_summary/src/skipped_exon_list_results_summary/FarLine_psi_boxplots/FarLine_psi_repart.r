#!/usr/bin/Rscript

psi_cols_select_fc <- function ( xxx ) {
  psi_cols <- colnames( xxx )[ grep( pattern='.*_psi$', x=colnames( xxx ) ) ]
  psi_cols <- c( 'id_gene', 'exon_skipped', psi_cols )

  return( xxx[ , psi_cols ] )
}

FarLine_psi_repart <- function( list_complete_tab ) {
  list_psi_tab <- lapply( list_complete_tab, psi_cols_select_fc )
  for ( xxx in names( list_psi_tab ) ) {
    write.table( list_psi_tab[[ xxx ]][ , c( 'id_gene', 'exon_skipped' ) ],
    file=paste0( "/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/FarLine_exons_results_summary/ref_exon_exprim/", xxx, '_exprim_list.tsv_temp' ),
    sep='_', col.names=FALSE, row.names=FALSE, quote=FALSE )
  }
}




# ##################
# #!/usr/bin/Rscript
#
# require( 'ggplot2' )
# require( 'reshape2' )
#
# intron_length_repart_plots_dir <- dirname(sys.frame(1)$ofile)
# source( paste( intron_length_repart_plots_dir, 'fig_functions_source.r', sep='/' ) )
# source( paste( intron_length_repart_plots_dir, 'fig_format.r', sep='/' ) )
# fig_functions_source( intron_length_repart_plots_dir )
#
# intron_length_repart_plots <- function ( final_data_list, nom_vec=c( 'query' ), out_dir, sign, color_pallette=NULL ) {
#   dir.create( out_dir, showWarnings=FALSE, recursive=TRUE )
#   col <- c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' )
#   if ( all( col %in% colnames( final_data_list[[ 1 ]] ) ) ) {
#
#     for ( col_name in col ) {
#       data_tab <- ggplot_repart_tab( final_data_list, col=col_name, nom_vec=nom_vec )
#       leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2
#
#       fig <- ggplot( data=data_tab, aes( x=value, y=repart, color=liste ) ) + scale_x_continuous( trans='log10' ) + ylab( 'fraction') + xlab( 'value' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
#       #~ coord_cartesian( xlim=c( 0, 10000 ) )
#
#       # set colors if defined
#       if ( ! is.null( color_pallette ) ) {
#         fig <- fig + scale_colour_manual( values=color_pallette )
#       }
#
#       png( fig_name( col_name, out_dir, sign, 'repart'), width=fig_width( col ) )
#       print(
#         fig + geom_point() + geom_line()
#       )
#       bouh <- dev.off()
#     }
#
#   }
#
#   # 2D density
#   data_tab <- lapply( final_data_list, function( xxx, col ) { xxx[ , col ] }, col=col )
#   data_tab <- do.call( rbind, data_tab )
#   data_tab$liste <- factor( rep( names( final_data_list ), unlist( lapply( final_data_list, nrow ) ) ), levels=nom_vec )
#
#   fig <- ggplot( data=data_tab, aes( x=INTRON_LENGTH_BEFORE, y=INTRON_LENGTH_AFTER, color=liste ) ) + scale_x_continuous( trans='log10' ) + scale_y_continuous( trans='log10' ) + ylab( 'INTRON_LENGTH_AFTER') + xlab( 'INTRON_LENGTH_BEFORE' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
#   #~ coord_cartesian( xlim=c( 0, 10000 ) )
#
#   # set colors if defined
#   if ( ! is.null( color_pallette ) ) {
#     fig <- fig + scale_colour_manual( values=color_pallette )
#   }
#
#   png( fig_name( 'intron_length', out_dir, sign, 'density2d'), width=fig_width( col ) )
#   print(
#     fig + geom_density2d()
#   )
#   bouh <- dev.off()
# }
