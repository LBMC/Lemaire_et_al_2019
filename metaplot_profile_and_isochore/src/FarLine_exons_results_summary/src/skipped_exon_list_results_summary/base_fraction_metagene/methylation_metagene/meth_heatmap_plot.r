#!/usr/bin/Rscript

meth_heatmap_plot_dir <- dirname(sys.frame(1)$ofile)
source( paste( meth_heatmap_plot_dir, 'base_fraction_metagene_fig_name.r', sep='/' ) )
source( paste( meth_heatmap_plot_dir, 'base_colors.r', sep='/' ) )
source( paste( meth_heatmap_plot_dir, 'sub_col_pal.r', sep='/' ) )

## build the curves of metagenes
meth_heatmap_plot_fc <- function( list_fraction_mat, prefix, vec_annot_name, suffix=NULL, scale_range=c( 0, 1 ) ) {
  # compute the name of the figure (and create the output folder)
  metagene_fig_name <- paste0( prefix, '_mCpG.png' )
  metagene_dir <- dirname( metagene_fig_name )
  heatmap_dir <- paste0( metagene_dir, '/heatmaps' )
  dir.create( heatmap_dir, showWarnings = FALSE, recursive = TRUE )

  for ( annot_name in vec_annot_name ) {
    out_prefix <- paste0( heatmap_dir, '/', unlist( strsplit( basename( metagene_fig_name ), split='.', fixed=TRUE ) )[ 1 ], '_', annot_name )
    if ( ! is.null( suffix ) ) {
      out_prefix <- paste0( out_prefix, '_', suffix )
    }
    fig_name <- paste0( out_prefix, '.png' )
    # message( fig_name )

    # compute the sum of signal
    ratio_vec <- apply( list_fraction_mat[[ annot_name ]], 1, function( xxx ) { sum( xxx, na.rm=TRUE )/( sum( !is.na( xxx ) ) ) } )

    # set the color palette
    col=base.colors( n=256, start=c( 0, 0, 1 ), end=c( 1, 0.65, 0 ) )
    col=sub_col_pal_fc( list_fraction_mat[[ annot_name ]], scale_range=scale_range, color_pallette=col )
    if ( NA %in% list_fraction_mat[[ annot_name ]] ) {
      list_fraction_mat[[ annot_name ]][ is.na( list_fraction_mat[[ annot_name ]] ) ] <- max( list_fraction_mat[[ annot_name ]], na.rm=TRUE ) + 1
      # col = c( col, '#000000FF' )
      col = c( col, '#FFFFFFFF' )
    }

    # build the heatmap
    png( fig_name, width=1200, height=1200 )
    heatmap( list_fraction_mat[[ annot_name ]][ rev( order( ratio_vec ) ), ], Rowv=NA, Colv=NA, scale='none', col=col, na.rm=FALSE )
    bouh <- dev.off()
  }
}
