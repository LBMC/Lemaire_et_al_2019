#!/usr/bin/Rscript

kmeans_heatmap <- function ( comp_pair, repl, suffix, out_dir, nbr_cluster, k_res_list, cov_df, cov_df_nullSig_list ) {
  # filename of the image
  # message( 'image' )
  if ( length( comp_pair ) == 1 ) {
    #~ prefix <- paste( c( comp_pair, ''), collapse='-' )
    sample <- paste( comp_pair, 'n', repl, '_', sep='' )
  } else if ( length( comp_pair ) == 2 ) {
    #~ prefix <- 'diffCov'
    sample <- paste( comp_pair[ 1 ], 'n', repl, '-', comp_pair[ 2 ], 'n', repl, '_', sep='' )
  }
  half_name <- paste( suffix, '.png' , sep='')
  heatmap_name <- paste( out_dir, '/heatmap_kmeansp2_', sample, 'allExons_', half_name, sep='' )

  # building the image
  couleurs <- rainbow( nbr_cluster + 1 )
  # cov_df_sub <- cov_df #~[ ! is.na( common_clusters ), ]
  common_clusters_sub <- k_res_list[[ repl ]]$cluster
  # cov_df_sub <- rbind( cov_df_sub[ order( common_clusters_sub ), ], cov_df_nullSig_list[[ repl ]] )
  cov_df <- rbind( cov_df[ order( common_clusters_sub ), ], cov_df_nullSig_list[[ repl ]] )
  common_clusters_sub <- c( common_clusters_sub, rep( nbr_cluster + 1, nrow( cov_df_nullSig_list[[ repl ]] ) ) )

  png( heatmap_name, width=1920, height=1200 )
  heatmap( as.matrix( cov_df ), Rowv=NA, RowSideColors=couleurs[ common_clusters_sub[ order( common_clusters_sub ) ] ], Colv=NA, scale='none', col=heat.colors( 256 ) )
  bouh <- dev.off()
}
