#!/usr/bin/Rscript

methTab2methMat_dir <- dirname(sys.frame(1)$ofile)
source( paste( methTab2methMat_dir, 'methTab2methVec.r', sep='/' ) )


## transform table of methylation for a set of annotation in a matrix of methylation rate per position and per annotation
methTab2methMat <- function( meth_tab, seq_tab, min.reads=1, start_pos=1, out_type='mat' ) {
  meth_mat <- list()

  for ( exon in 1:nrow( seq_tab ) ) {
    seq_line <- seq_tab[ exon, ]

    # parse the coordinates
    coord_parts <- strsplit( seq_line[ 1, 'coordinates' ], split=':' )[[ 1 ]]
    chrom <- coord_parts[ 1 ]
    borders <- coord_parts[ 2 ]

    borders_parts <- strsplit( borders, split='-' )[[ 1 ]]
    start <- as.numeric( borders_parts[ 1 ] )
    end <- as.numeric( borders_parts[ 2 ] )

    # build vector of methylation
    meth_vec <- methTab2methVec( meth_tab, chrom, start, end, min.reads=min.reads )
    if ( seq_line[ 1, 'strand' ] == -1 ) {
      meth_vec <- rev( meth_vec )
    }
    names( meth_vec ) <- as.character( start_pos:( start_pos + length( meth_vec ) - 1 ) )

    meth_mat[[ exon ]] <- meth_vec
  }

  if ( out_type == 'mat' ) {
    meth_mat <- do.call( rbind, meth_mat )
  }

  return( meth_mat )
}
