#!/usr/bin/Rscript

###########
listMeta2DF <- function( list_metagene ) {
  # reshape the list in a data frame
  df_metagene <- melt( list_metagene )

  # add the position for each base
  col_name <- 'position'
  df_metagene[ col_name ] <- 0 #~ rep( 0, nrow( df_metagene ) )
  for ( xxx in names( list_metagene ) ) {
    df_metagene[ df_metagene$L1 == xxx, col_name ] <- as.numeric( names( list_metagene[[ xxx ]] ) )
  }

  return( df_metagene )
}
