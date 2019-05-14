#!/usr/bin/Rscript

require( 'bigWig' )
stranded_coverages_varLen_dir <- paste( dirname(sys.frame(1)$ofile), 'metagene_coverage_functions', sep='/' )

## absolute coverage in one condition
stranded_coverages_varLen <- function( bw_file, annot, nbins=10, off_set=1, normFactor=1, harmonization = FALSE, ref_add_chr=F, binning_fc=mean, gene_concat=FALSE ) {
  # bw_reader <- load.bigWig( filename=bw_file )

  ## check if human chromosome
  if ( ref_add_chr ) {
    annot[ , 1 ] <- paste0( 'chr', annot[ , 1 ] )
  }

  # # retrieve coverage from bam file for all annotations and apply normalization
  # cov_df <- matrix( NA, nrow=nrow( annot ), ncol=nbins, dimnames=list( annot[ , 4 ], 1:nbins ) )
  #
  # cov_list <- list()
  # for ( annot_idx in 1:nrow( annot ) ) {
  #   cov_annot <- step.probeQuery.bigWig( bw=bw_reader, chrom=annot[ annot_idx, 1 ], start=annot[ annot_idx, 2 ], end=annot[ annot_idx, 3 ], step=1, gap.value=0 )
  #   if ( annot[ annot_idx, 6 ] == '-1' ) {
  #     cov_annot <- rev( cov_annot )
  #   }
  #   ends <- ceiling( ( ( 1:nbins ) / nbins ) * length( cov_annot ) )
  #   start <- c( 1, head( ends, n=-1 ) )
  #   for ( bin_idx in 1:length( start ) ) {
  #     cov_df[ annot_idx, bin_idx ] <- binning_fc( cov_annot[ start[ bin_idx ]:ends[ bin_idx ] ] )
  #   }
  # }

  ##
  script <- paste0( stranded_coverages_varLen_dir, '/retrieve_bw_values_binned.py' )
  # script <- paste0( stranded_coverages_varLen_dir, '/retrieve_bw_values_binned' )
  annot_file <- tempfile()
  write.table( annot, file=annot_file, quote=F, sep="\t", col.names=F, row.names=F)
  #~ write.table( annot, file='./test.bed', quote=F, sep="\t", col.names=F, row.names=F)

  out_file <- tempfile()

  cmd=paste( 'python3', script, bw_file, annot_file, out_file, nbins, sep=' ' )
  if ( gene_concat ) {
    cmd=paste( cmd, '--group-concat', sep=' ' )
  }
  # cmd=paste( script, bw_file, annot_file, out_file, nbins, sep=' ' )
  message( cmd )
  bouh <- system( cmd )
  cov_df <- as.matrix( read.table( out_file, h=F, stringsAsFactors=F ) )
  colnames( cov_df ) = 1:nbins

  ##


  if ( ! harmonization ) {
    cov_df <- cov_df / normFactor
  }

  start_pos=( 1 - off_set )
  colnames( cov_df ) <- seq( from=start_pos, to=( start_pos + ncol( cov_df ) - 1 )  )

  # harmonize the scale among the exons
  if ( harmonization ) {
    idx <- which( rowSums( cov_df ) != 0 )
    cov_df[ idx, ] <- as.data.frame( t( apply( cov_df[ idx, ], 1, function( xxx ) { xxx / max( xxx ) } ) ) )
  }

  return( cov_df )
}
