#!/usr/bin/Rscript

pattern_finder_fc_dir <- dirname(sys.frame(1)$ofile)
source( paste( pattern_finder_fc_dir, 'regexpr_finder_fc.r', sep='/' ) )

###########
pattern_finder_fc <- function ( sequence, base_vec, skew_measure=NULL, use.regexpr=FALSE, all.pos=TRUE ) {

  if ( use.regexpr ) {
    # use regular expression engine to find the patterns
    logic_seq <- rep( FALSE, length( sequence )  )
    for ( pattern in base_vec ) {
      bouh <- regexpr_finder_fc( sequence, pattern )
      logic_seq <- logic_seq | bouh
    }

  } else {
    # custom search engine to find the patterns

    # evaluate the length of the patterns to look for
    base_vec_len <- unlist( lapply( strsplit( base_vec, split='' ), length ) )
    logic_seq <- ( sequence %in% base_vec[ base_vec_len == 1 ] )

    # count the occurences, a way according to the length of the patterns to look for
    if ( any( base_vec_len != 1 ) ) {

      eff_seq <- sequence[ 1:( length( sequence ) - max( base_vec_len ) + 1 ) ]
      index_mat <- matrix( c( 1:length( eff_seq ) ), ncol=max( base_vec_len ), nrow=length( eff_seq ) )
      for ( xxx in c( 2:max( base_vec_len ) ) ) {
        index_mat[ , xxx ] <- index_mat[ , xxx - 1 ] + 1

        eff_seq <- paste( eff_seq, sequence[ xxx:( length( sequence ) - max( base_vec_len ) + xxx ) ], sep='' )
        eff_logic_seq <- eff_seq %in% base_vec[ base_vec_len == xxx ]

        col_idx <- 1
        if ( all.pos ) {
          col_idx <- c( 1:xxx )
        }
        index_logic_seq <- unique( as.vector( index_mat[ eff_logic_seq, col_idx ] ) )

        logic_seq[ index_logic_seq ] <- TRUE
      }
    }
  }

  return( logic_seq )
}
