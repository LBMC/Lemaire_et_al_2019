#!/usr/bin/Rscript

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)
list_exon <- read.table( args[ 1 ], h=T, stringsAsFactors=F, sep='\t' )
splicing_site <- args[ 2 ]
in_exon <- as.numeric( args[ 3 ] )
out_off <- as.numeric( args[ 4 ] )
exon <- read.table( args[ 5 ], h=T, stringsAsFactors=F, sep='\t' )

# recover lines of firts and last exons ( to avoid computing splicing sites on first/last exons )
exon[ 'egp' ] <- paste0( exon$id_gene, '_', exon$pos_sur_gene )
first_exons <- tapply( exon$pos_sur_gene, exon$id_gene, min )
last_exons <- tapply( exon$pos_sur_gene, exon$id_gene, max )
first_exons <- paste0( names( first_exons ), '_', first_exons )
last_exons <- paste0( names( last_exons ), '_', last_exons )
border_exons <- union( first_exons, last_exons )

# treat coordinates
list_exon[ 'chrom' ] <- apply( list_exon[ 'coord' ], 1, function( xxx ) { strsplit( xxx, split=':' )[[ 1 ]][ 1 ] } )
list_exon[ 'start' ] <- as.numeric( apply( list_exon[ 'coord' ], 1, function( xxx ) { strsplit( strsplit( xxx, split=':' )[[ 1 ]][ 2 ], split='-' )[[ 1 ]][ 1 ] } ) )
list_exon[ 'end' ] <- as.numeric( apply( list_exon[ 'coord' ], 1, function( xxx ) { strsplit( strsplit( xxx, split=':' )[[ 1 ]][ 2 ], split='-' )[[ 1 ]][ 2 ] } ) )

# compute borders of splicing site regions
ss_tab <- list_exon
ss_tab$coord <- ':-'
ss_tab$start <- -1
ss_tab$end <- -1
if ( splicing_site == '3SS' ) {
  pos_str <- ss_tab$strand == 1
  ss_tab$start[ pos_str ] <- list_exon$start[ pos_str ] - out_off
  ss_tab$end[ pos_str ] <- list_exon$start[ pos_str ] + in_exon - 1

  neg_str <- ss_tab$strand == -1
  ss_tab$start[ neg_str ] <- list_exon$end[ neg_str ] - in_exon + 1
  ss_tab$end[ neg_str ] <- list_exon$end[ neg_str ] + out_off

} else if ( splicing_site == '5SS' ) {
  pos_str <- ss_tab$strand == 1
  ss_tab$start[ pos_str ] <- list_exon$end[ pos_str ] - in_exon + 1
  ss_tab$end[ pos_str ] <- list_exon$end[ pos_str ] + out_off

  neg_str <- ss_tab$strand == -1
  ss_tab$start[ neg_str ] <- list_exon$start[ neg_str ] - out_off
  ss_tab$end[ neg_str ] <- list_exon$start[ neg_str ] + in_exon - 1

}

ss_tab$coord <- paste0( ss_tab$chrom, ':', ss_tab$start, '-', ss_tab$end )

# look for column and type of exon id
id_col <- 'exon_id'
if ( 'exon_gene_pos' %in% colnames( ss_tab ) ) {
  colnames( ss_tab )[ colnames( ss_tab ) == 'exon_gene_pos' ] <- 'exon_id'
}

# filter out the first and last exons
if ( ! '--coord-only' %in% args ) {
  if ( splicing_site == '3SS' ) {
    ss_tab <- ss_tab[ ! ss_tab[ , id_col ] %in% first_exons, ]
  } else if ( splicing_site == '5SS' ) {
    ss_tab <- ss_tab[ ! ss_tab[ , id_col ] %in% last_exons, ]
  }
}


# output in stdout the table of splicing site regions
out_col <- c( id_col, 'coord', 'strand' )
write.table( ss_tab[ out_col ], file=stdout(), col.names=T, row.names=F, sep='\t', quote=F )
