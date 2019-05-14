#!/usr/bin/Rscript

require( 'bigWig' )

#### retrieve sequencing depth for one samples
bw_seqDepth_retriever <- function ( bw_file ) {
  bw_reader <- load.bigWig( bw_file, udcDir='./udcCache' )
  return( bw_reader$basesCovered * bw_reader$mean )
}
