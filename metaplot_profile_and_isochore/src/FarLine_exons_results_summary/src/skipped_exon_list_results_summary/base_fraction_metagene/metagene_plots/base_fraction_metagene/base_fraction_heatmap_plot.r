#!/usr/bin/Rscript

base_fraction_heatmap_plot_dir <- dirname(sys.frame(1)$ofile)
source( paste( base_fraction_heatmap_plot_dir, 'base_fraction_metagene_fig_name.r', sep='/' ) )
source( paste( base_fraction_heatmap_plot_dir, 'base_colors.r', sep='/' ) )
source( paste( base_fraction_heatmap_plot_dir, 'sub_col_pal.r', sep='/' ) )

## build the curves of metagenes
base_fraction_heatmap_plot_fc <- function( list_fraction_mat, vec_annot_name, base_vec, color_pallette=NULL, scale_range=c( 0, 1 ), skew_measure=NULL, prefix='', entropy=FALSE, struct_info=FALSE, ent_base_vec=NULL, regexpr_fig_name=NULL, presence.only=FALSE ) {
  for ( annot_name in vec_annot_name ) {
    # compute the name of the figure (and create the output folder)
    metagene_fig_name <- base_fraction_metagene_fig_name_fc( base_vec, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )
    metagene_dir <- dirname( metagene_fig_name )
    heatmap_dir <- paste0( metagene_dir, '/heatmaps' )
    dir.create( heatmap_dir, showWarnings = FALSE, recursive = TRUE )

    out_prefix <- paste0( prefix, '_', annot_name )
    out_prefix <- paste0( heatmap_dir, '/', unlist( strsplit( basename( metagene_fig_name ), split='.', fixed=TRUE ) )[ 1 ], '_', annot_name )
    fig_name <- paste0( out_prefix, '.png' )
    # message( fig_name )

    # compute the sum of signal
    sum_vec <- apply( list_fraction_mat[[ annot_name ]], 1, sum )

    # set the color palette
    col=base.colors( n=256, start=c( 0, 0, 1 ), end=c( 1, 0.65, 0 ) )
    col=sub_col_pal_fc( list_fraction_mat[[ annot_name ]], scale_range=scale_range, color_pallette=col )
    if ( NA %in% list_fraction_mat[[ annot_name ]] ) {
      list_fraction_mat[[ annot_name ]][ is.na( list_fraction_mat[[ annot_name ]] ) ] <- max( list_fraction_mat[[ annot_name ]], na.rm=TRUE ) + 1
      col = c( col, '#FFFFFFFF' )
    }

    # build the heatmap
    png( fig_name, width=1200, height=1200 )
    heatmap( list_fraction_mat[[ annot_name ]][ rev( order( sum_vec ) ), ], Rowv=NA, Colv=NA, scale='none', col=col, na.rm=FALSE )
    bouh <- dev.off()
  }
}
