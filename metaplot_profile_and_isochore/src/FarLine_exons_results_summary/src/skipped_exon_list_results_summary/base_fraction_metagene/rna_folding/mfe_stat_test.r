#!/usr/bin/Rscript

mfe_stat_test <- function( mfe_df_list, names_vec, stat_test_fc=wilcox.test, ... ) {
  write( 'Test difference in MFE', file=stdout() )
  if ( length( names_vec ) == 1 ) {
    write( '!!! No statistical test to do. Only one liste.', file=stdout() )
  } else {
    for ( name_idx in 1:( length( names_vec ) - 1 ) ) {
      for ( name_idx_2 in ( name_idx + 1 ):length( names_vec ) ) {
        write( paste0( names_vec[ name_idx ], ' _vs_ ', names_vec[ name_idx_2 ] ), file=stdout() )
        p_val <- stat_test_fc( mfe_df_list[[ names_vec[ name_idx ] ]][[ 'MFE' ]], mfe_df_list[[ names_vec[ name_idx_2 ] ]][[ 'MFE' ]], ... )[[ 'p.value' ]]
        write( paste0( 'p-value: ', p_val ), file=stdout() )
      }
    }
  }

  ####
}
