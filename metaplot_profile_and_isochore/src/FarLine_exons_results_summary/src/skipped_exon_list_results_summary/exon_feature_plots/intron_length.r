#!/usr/bin/Rscript

intron_length_dir <- dirname(sys.frame(1)$ofile)
source( paste( intron_length_dir, 'val_retriever.r', sep='/' ) )
source( paste( intron_length_dir, 'flanquing_exon_feat.r', sep='/' ) )

# compute the intron length (up or down of the skipped exon)
intron_length <- function ( exon_tab_complete, feat_tab, pos ) {
  # exon_tab_complete needs the columns: 'START.END', 'exons_flanquants', 'id_gene'
  # feat_tab needs the columns: 'START.END', 'join_id'
  if ( pos == 'up' ) { fields <- c( 1, 2 ) } else if ( pos == 'down' ) { fields <- c( 2, 1 ) }

  exon_border <- as.numeric( val_retriever( exon_tab_complete[ , 'START.END' ], sep='-', field=fields[ 1 ] ) )
  voisin <- paste( exon_tab_complete$id_gene, val_retriever( exon_tab_complete[ , 'exons_flanquants' ], sep=':', field=fields[ 1 ] ), sep='_' )
  voisin_border <- flanquing_exon_feat( exon_tab_complete, feat_tab, pos=pos )[ , 'START.END' ]
  voisin_border <- as.numeric( val_retriever( voisin_border, sep='-', field=fields[ 2 ] ) )
  abs( exon_border - voisin_border ) - 1

}
