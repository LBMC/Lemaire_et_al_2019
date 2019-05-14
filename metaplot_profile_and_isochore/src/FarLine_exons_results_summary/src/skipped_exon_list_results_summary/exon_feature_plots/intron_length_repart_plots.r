#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

intron_length_repart_plots_dir <- dirname(sys.frame(1)$ofile)
source( paste( intron_length_repart_plots_dir, 'fig_functions_source.r', sep='/' ) )
source( paste( intron_length_repart_plots_dir, 'fig_format.r', sep='/' ) )
fig_functions_source( intron_length_repart_plots_dir )

intron_length_repart_plots <- function ( final_data_list, nom_vec=c( 'query' ), out_dir, sign, color_pallette=NULL ) {
  ## add columns for the min and max flanquing intron length
  final_data_list <- lapply( final_data_list, function ( xxx ) {
    xxx[ 'INTRON_LENGTH_MIN' ] <- apply( xxx[ c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' ) ], 1, min, na.rm=TRUE )
    xxx[ 'INTRON_LENGTH_MAX' ] <- apply( xxx[ c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' ) ], 1, max, na.rm=TRUE )
    return( xxx )
  } )

  dir.create( out_dir, showWarnings=FALSE, recursive=TRUE )
  col <- c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER', 'INTRON_LENGTH_MIN', 'INTRON_LENGTH_MAX' )
  if ( all( col %in% colnames( final_data_list[[ 1 ]] ) ) ) {

    for ( col_name in col ) {
      data_tab <- ggplot_repart_tab( final_data_list, col=col_name, nom_vec=nom_vec )
      leg_nrow <- ( length( levels( factor( data_tab$liste ) ) ) + 1 ) %/% 2

      fig <- ggplot( data=data_tab, aes( x=value, y=repart, color=liste ) ) + scale_x_continuous( trans='log10' ) + ylab( 'fraction') + xlab( 'value' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
      #~ coord_cartesian( xlim=c( 0, 10000 ) )

      # set colors if defined
      if ( ! is.null( color_pallette ) ) {
        fig <- fig + scale_colour_manual( values=color_pallette )
      }

      png( fig_name( col_name, out_dir, sign, 'repart') )
      print(
        fig + geom_point() + geom_line()
      )
      bouh <- dev.off()
    }

  }

  # 2D density of exons for their flanquing intron lengths (upstream and downstream)
  data_tab <- lapply( final_data_list, function( xxx, col ) { xxx[ , col ] }, col=col )
  data_tab <- do.call( rbind, data_tab )
  data_tab$liste <- factor( rep( names( final_data_list ), unlist( lapply( final_data_list, nrow ) ) ), levels=nom_vec )

  fig <- ggplot( data=data_tab, aes( x=INTRON_LENGTH_BEFORE, y=INTRON_LENGTH_AFTER, color=liste ) ) + scale_x_continuous( trans='log10' ) + scale_y_continuous( trans='log10' ) + ylab( 'INTRON_LENGTH_AFTER') + xlab( 'INTRON_LENGTH_BEFORE' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
  #~ coord_cartesian( xlim=c( 0, 10000 ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    fig <- fig + scale_colour_manual( values=color_pallette )
  }

  png( fig_name( 'intron_length', out_dir, sign, 'density2d') )
  print(
    fig + geom_density2d()
  )
  bouh <- dev.off()

  # 2D density of exons for their flanquing intron lengths (min and max)
  data_tab <- lapply( final_data_list, function( xxx, col ) { xxx[ , col ] }, col=col )
  data_tab <- do.call( rbind, data_tab )
  data_tab$liste <- factor( rep( names( final_data_list ), unlist( lapply( final_data_list, nrow ) ) ), levels=nom_vec )

  fig <- ggplot( data=data_tab, aes( x=INTRON_LENGTH_MIN, y=INTRON_LENGTH_MAX, color=liste ) ) + scale_x_continuous( trans='log10' ) + scale_y_continuous( trans='log10' ) + ylab( 'INTRON_LENGTH_MAX') + xlab( 'INTRON_LENGTH_MIN' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
  #~ coord_cartesian( xlim=c( 0, 10000 ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    fig <- fig + scale_colour_manual( values=color_pallette )
  }

  png( fig_name( 'intron_length_min-max', out_dir, sign, 'density2d') )
  print(
    fig + geom_density2d()
  )
  bouh <- dev.off()


  # #### Filter for intron under 2000
  # col <- c( 'join_id', 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' )
  # sbh_fc <- function() { geom_hex( bins=60 ) }
  #
  # ####
  # png( fig_name( 'intron_length_dnmt3b', out_dir, sign, 'density2d' ), width=fig_width( col ) )
  # fig_3b <- ggplot( data=data_tab[ data_tab[ , 'liste' ] == 'exon_down_siDNMT3b-siGL2', ], aes( x=INTRON_LENGTH_BEFORE, y=INTRON_LENGTH_AFTER ) ) + scale_x_continuous( trans='log10' ) + scale_y_continuous( trans='log10' ) + ylab( 'INTRON_LENGTH_AFTER') + xlab( 'INTRON_LENGTH_BEFORE' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
  # print(
  #   fig_3b + sbh_fc()
  # )
  # bouh <- dev.off()
  #
  # small_intron <- data_tab[ data_tab[ , 'liste' ] == 'exon_down_siDNMT3b-siGL2' & data_tab[ , 'INTRON_LENGTH_BEFORE' ] <= 2000 & data_tab[ , 'INTRON_LENGTH_AFTER' ] <= 2000, ]
  # print( str( small_intron ) )
  # write.table( small_intron[ , 'join_id' ], file='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/FarLine_exons_results_summary/siDNMT3b-siGL2/subGroups/exon_down_siDNMT3b-siGL2_smallIntron.tsv', quote=FALSE, col.names=TRUE, row.names=FALSE )
  #
  # max_intron <- apply( data_tab[ data_tab[ , 'liste' ] == 'exon_down_siDNMT3b-siGL2', c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' ) ], 1, max )
  # png( fig_name( 'intron_length_dnmt3b', out_dir, sign, 'repartMaxIntron' ), width=fig_width( col ) )
  # plot( log2( max_intron[ order( max_intron ) ] ), c( 1:length( max_intron ) )/length( max_intron ), pch=16 )
  # grid()
  # bouh <- dev.off()
  #
  # max_intron_filt <- apply( small_intron[ , c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' ) ], 1, max )
  # png( fig_name( 'intron_length_dnmt3b_mI2000', out_dir, sign, 'repartMaxIntron' ), width=fig_width( col ) )
  # plot( log2( max_intron_filt[ order( max_intron_filt ) ] ), head( c( 1:length( max_intron ) )/length( max_intron ), n=length( max_intron_filt ) ), pch=16 )
  # grid()
  # bouh <- dev.off()
  #
  #
  # ####
  # png( fig_name( 'intron_length_GC_55', out_dir, sign, 'density2d' ), width=fig_width( col ) )
  # fig_3b <- ggplot( data=data_tab[ data_tab[ , 'liste' ] == 'GC_55-65', ], aes( x=INTRON_LENGTH_BEFORE, y=INTRON_LENGTH_AFTER ) ) + scale_x_continuous( trans='log10' ) + scale_y_continuous( trans='log10' ) + ylab( 'INTRON_LENGTH_AFTER') + xlab( 'INTRON_LENGTH_BEFORE' ) + fig_format + guides(fill = guide_legend(nrow = leg_nrow))
  # print(
  #   fig_3b + sbh_fc()
  # )
  # bouh <- dev.off()
  #
  # small_intron <- data_tab[ data_tab[ , 'liste' ] == 'GC_55-65' & data_tab[ , 'INTRON_LENGTH_BEFORE' ] <= 2000 & data_tab[ , 'INTRON_LENGTH_AFTER' ] <= 2000, ]
  # print( str( small_intron ) )
  #
  # max_intron <- apply( data_tab[ data_tab[ , 'liste' ] == 'GC_55-65', c( 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER' ) ], 1, max )
  # png( fig_name( 'intron_length_GC_55', out_dir, sign, 'repartMaxIntron' ), width=fig_width( col ) )
  # plot( log2( max_intron[ order( max_intron ) ] ), c( 1:length( max_intron ) )/length( max_intron ), pch=16 )
  # grid()
  # bouh <- dev.off()

}
