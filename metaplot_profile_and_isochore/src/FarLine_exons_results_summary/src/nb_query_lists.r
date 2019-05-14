#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
exon_tables <- strsplit( args[ 1 ], split='+', fixed=TRUE )[[ 1 ]]
exon_tables_nb <- strsplit( args[ 2 ], split='+', fixed=TRUE )[[ 1 ]]
nb_type <- args[ 3 ]
exon_fdb <- args[ 4 ]

# set modification on exon id/pos according to neighbour
if ( nb_type == 'nbup' ) {
  trans_val <- -1
  fl_col <- 1
} else if ( nb_type == 'nbdown' ) {
  trans_val <- 1
  fl_col <- 2
}

# load the exons_genomiques_bis table
exon_fdb_tab <- read.table( exon_fdb, header=TRUE, stringsAsFactors=FALSE )
exon_fdb_tab[ 'exon_gene_pos' ] <- paste( exon_fdb_tab$id_gene, exon_fdb_tab$pos_sur_gene, sep='_' )
exon_fdb_tab[ 'id_exon' ] <- exon_fdb_tab[ 'id' ]
## fix unordered coordinates
coord_fdb <- exon_fdb_tab[ c( 'start_sur_chromosome', 'end_sur_chromosome' ) ]
exon_fdb_tab[ 'start_sur_chromosome' ] <- apply( coord_fdb, 1, min )
exon_fdb_tab[ 'end_sur_chromosome' ] <- apply( coord_fdb, 1, max )

exon_fdb_tab[ 'coordonnees' ] <- paste0( exon_fdb_tab[ , 'chromosome' ], ':', exon_fdb_tab[ , 'start_sur_chromosome' ], '-', exon_fdb_tab[ , 'end_sur_chromosome' ] )

for ( idx in 1:length( exon_tables ) ) {
  exon_tab <- exon_tables[ idx ]

  # load the all_data table
  adt_nb <- read.table( exon_tab, header=TRUE, stringsAsFactors=FALSE )
  out_col <- colnames( adt_nb )
  exon_pos_query <- adt_nb[ , 'exon_pos' ]

  ## update the exon position
  exon_fl <- do.call( rbind, strsplit( adt_nb$exons_flanquants, split=':' ) )
  adt_nb[ 'exon_pos' ] <- as.numeric( exon_fl[ , fl_col ] )

  ## update the flanquing exons
  exon_fl[ , 3 - fl_col ] <- exon_pos_query
  exon_fl[ , fl_col ] <- as.character( as.numeric( exon_fl[ , fl_col ] ) + trans_val )
  adt_nb[ 'exons_flanquants' ] <- paste( exon_fl[ , 1 ], exon_fl[ , 2 ], sep=':' )

  ## recover the coordinates of the neighbour exon, and the exon ids
  adt_nb[ 'exon_gene_pos' ] <- paste( adt_nb$id_gene, adt_nb$exon_pos, sep='_' )
  adt_nb <- merge( adt_nb, exon_fdb_tab[ c( 'exon_gene_pos', 'id_exon', 'coordonnees', 'strand' ) ], by='exon_gene_pos', suffixes=c( '.query', '' ), sort=FALSE )

  ## write the transformed table
  dir.create( dirname( exon_tables_nb[ idx ] ), recursive=TRUE, showWarnings=FALSE )
  write.table( adt_nb[ out_col ], exon_tables_nb[ idx ], quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t' )
}
