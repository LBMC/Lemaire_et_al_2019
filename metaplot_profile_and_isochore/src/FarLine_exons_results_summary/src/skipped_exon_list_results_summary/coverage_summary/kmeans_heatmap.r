#!/usr/bin/Rscript

kmeans_heatmap <- function( heatmap_name, nbr_cluster, k_res, cov_df, cov_df_nullSig, ymax=NULL ) {
  couleurs <- rainbow( nbr_cluster + 1 )
  # cov_df_sub <- cov_df #~[ ! is.na( common_clusters ), ]
  common_clusters_sub <- k_res_list[[ repl ]]$cluster
  # cov_df_sub <- rbind( cov_df_sub[ order( common_clusters_sub ), ], cov_df_nullSig_list[[ repl ]] )
  cov_df <- rbind( cov_df[ order( common_clusters_sub ), ], cov_df_nullSig_list[[ repl ]] )
  common_clusters_sub <- c( common_clusters_sub, rep( nbr_cluster + 1, nrow( cov_df_nullSig_list[[ repl ]] ) ) )

  cov_df <- as.matrix( cov_df )
  if ( ! is.null( ymax ) ) {
    cov_df[ cov_df > ymax ] <- ymax
  }

  png( heatmap_name, width=1200, height=1200 )
  heatmap( cov_df, Rowv=NA, RowSideColors=couleurs[ common_clusters_sub[ order( common_clusters_sub ) ] ], Colv=NA, scale='none', col=heat.colors( 256 ) )
  dev.off()
}
