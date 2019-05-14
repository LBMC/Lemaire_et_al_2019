#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

feat_fig_builder_dir <- dirname(sys.frame(1)$ofile)
source( paste( feat_fig_builder_dir, 'intron_length_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'intron_length_repart_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'exon_length_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'exon_length_repart_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'acceptor_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'donor_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'diffAcceptor_plots.r', sep='/' ) )
source( paste( feat_fig_builder_dir, 'diffDonor_plots.r', sep='/' ) )


feat_fig_builder <- function ( final_data_list, nom_vec, nom_query, out_dir='./', sign, color_pallette=NULL ) {

  dir.create( out_dir, showWarnings=FALSE )
  draw_quantiles=c( 0.25, 0.5, 0.75 )

  ##
  # intron length
  intron_length_plots( final_data_list, nom_vec, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )
  intron_length_repart_plots( final_data_list, nom_vec, paste0( out_dir, '/reparts' ), sign, color_pallette=color_pallette )

  ##
  # exon_length
  exon_length_plots( final_data_list, nom_vec, nom_query, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )
  exon_length_repart_plots( final_data_list, nom_vec, nom_query, paste0( out_dir, '/reparts' ), sign, color_pallette=color_pallette )

  ##
  # acceptor (3'SS) force
  acceptor_plots( final_data_list, nom_vec, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )

  ##
  # donor (5'SS) force
  donor_plots( final_data_list, nom_vec, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )

  ##
  # acceptor (3'SS) force
  diffAcceptor_plots( final_data_list, nom_vec, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )

  ##
  # donor (5'SS) force
  diffDonor_plots( final_data_list, nom_vec, out_dir, sign, draw_quantiles=draw_quantiles, color_pallette=color_pallette )


}
