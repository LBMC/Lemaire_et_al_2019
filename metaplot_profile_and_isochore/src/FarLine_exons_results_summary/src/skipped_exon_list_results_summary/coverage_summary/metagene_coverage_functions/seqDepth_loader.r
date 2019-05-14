#!/usr/bin/Rscript

#### loader of seqDepth file
seqDepth_loader <- function( seqDepth_file ) {
  seqDepth_df <- read.table( seqDepth_file, stringsAsFactor=FALSE, header=FALSE )
  colnames( seqDepth_df ) <- c( 'condition', 'replicate', 'sequencing_depth' )
  seqDepth_df <- data.frame( seqDepth_df, normFactor=seqDepth_df[ ,'sequencing_depth' ]/min( seqDepth_df[ ,'sequencing_depth' ] ) )
  return( seqDepth_df )
}
