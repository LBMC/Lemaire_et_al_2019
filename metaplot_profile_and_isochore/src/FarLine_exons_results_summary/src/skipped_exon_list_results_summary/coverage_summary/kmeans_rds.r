#!/usr/bin/Rscript

kmeans_rds <- function ( out_RDS_dir, suffix, k_res_list, cov_enough_id, min_max_signal, annot_file, nbr_cluster, comp_pair ) {
  if ( length( comp_pair ) == 1 ) {
    prefix <- comp_pair
  } else if ( length( comp_pair ) == 2 ) {
    prefix <- paste( comp_pair, collapse='-' )
  }
  saveRDS( object=list(
  k_res_list=k_res_list,
  cov_enough_id=cov_enough_id,
  min_max_signal=min_max_signal,
  annot_file=annot_file,
  nbr_cluster=nbr_cluster,
  condition=comp_pair
  ), file=paste( out_RDS_dir, '/cluster_data_', prefix, '_', suffix, '.RDS', sep='' ) )
}
