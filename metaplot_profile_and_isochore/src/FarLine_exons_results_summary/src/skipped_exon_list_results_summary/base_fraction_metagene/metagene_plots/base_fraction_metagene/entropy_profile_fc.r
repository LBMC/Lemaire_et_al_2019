#!/usr/bin/Rscript

entropy_profile_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( entropy_profile_fc_dir, 'entropy_fc.r', sep='/' ) )


###########
entropy_profile_fc <- function ( counts_mat, struct_info=FALSE ) {
  ##> struct_info, when TRUE, specify to retrieve the entropy to the theoritical maximum entropy, giving a measure of "structuration, information" from the counts

  entropy_profile <- apply( counts_mat, 2, entropy_fc )
  if ( struct_info ) {
    max_ent <- entropy_fc( rep( 1, nrow( counts_mat ) ) )
    entropy_profile <- max_ent - entropy_profile
  }

  return( entropy_profile )
}
