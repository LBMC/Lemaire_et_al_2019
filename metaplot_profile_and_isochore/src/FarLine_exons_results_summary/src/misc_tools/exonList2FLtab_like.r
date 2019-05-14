#!/usr/bin/Rscript

all_args <- commandArgs(trailingOnly = F)
script_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

exon_FasterDB_tab <- paste0( script_dir, '/../../data/exons_genomiques_bis.tsv' )
exon <- read.table( exon_FasterDB_tab, h=T, stringsAsFactors=F )
exon[ 'exon_gene_pos' ] <- paste0( exon$id_gene, '_', exon$pos_sur_gene )
exon[ 'feat_coord' ] <- paste0(
  exon$chromosome,
  ':',
  apply( exon[ c( 'start_sur_chromosome', 'end_sur_chromosome' ) ], 1, min ),
  '-',
  apply( exon[ c( 'start_sur_chromosome', 'end_sur_chromosome' ) ], 1, max )
)

args <- commandArgs(trailingOnly = TRUE)
exon_list_path <- args[ 1 ]
header=T
if ( '--no-header' %in% args ) {
  header=F
}
exon_query <- read.table( exon_list_path, h=header, stringsAsFactors=F )
if ( ! header ) {
  colnames( exon_query ) <- 'exon_gene_pos'
}
exon_query[ 'id_event' ] <- 1:nrow( exon_query )

gene_path <- paste0( script_dir, '/../../data/genes.tsv' )
gene <- read.table( gene_path, h=T, stringsAsFactors=F )

# print( str( exon_query ), file=stderr() )
temp_mat <- do.call( rbind, strsplit( exon_query$exon_gene_pos, split='_' ) )
exon_query[ 'exon_pos' ] <- as.numeric( temp_mat[ , 2 ] )
exon_query[ 'id_gene' ] <- as.numeric( temp_mat[ , 1 ] )
exon_query[ 'exons_flanquants' ] <- paste0( exon_query$exon_pos - 1, ':', exon_query$exon_pos + 1 )

exon_query <- merge( exon_query, gene[ c( 'id', 'official_symbol' ) ], by.x='id_gene', by.y='id', sort=F, all.x=T )
exon_query[ 'gene_symbol' ] <- exon_query[ 'official_symbol' ]

exon_query <- merge( exon_query, exon[ c( 'exon_gene_pos', 'id' ) ], sort=F, all.x=T )
exon_query[ 'id_exon' ] <- exon_query[ 'id' ]

exon_query <- merge( exon_query, exon[ c( 'exon_gene_pos', 'feat_coord' ) ], sort=F, all.x=T )
exon_query[ 'coordonnees' ] <- exon_query[ 'feat_coord' ]

exon_query <- exon_query[ order( exon_query$id_event ), ]

out_col <- c( 'id_event', 'id_exon', 'exons_flanquants', 'coordonnees', 'id_gene', 'gene_symbol', 'exon_pos' )
write.table( exon_query[ out_col ], file=stdout(), quote=F, sep='\t', row.names=F, col.names=T )

## id_event
## id_exon
## exons_flanquants
## coordonnees
## id_gene
## gene_symbol
## exon_pos


####
