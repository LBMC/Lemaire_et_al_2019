#!/usr/bin/Rscript

######
pos_SS <- function ( bed_tab, splicingSite ) {
  bed_tab[ , 'end' ] <- bed_tab[ , 'end' ] + 1
  SS_vec <- rep( -1, nrow( merge_tab ) )

  if ( splicingSite == '3SS' ) {
    SS_vec[ bed_tab$strand == 1 ] <- bed_tab$start[ bed_tab$strand == 1 ]
    SS_vec[ bed_tab$strand == -1 ] <- bed_tab$end[ bed_tab$strand == -1 ]
  } else if ( splicingSite == '5SS' ) {
    SS_vec[ bed_tab$strand == 1 ] <- bed_tab$end[ bed_tab$strand == 1 ]
    SS_vec[ bed_tab$strand == -1 ] <- bed_tab$start[ bed_tab$strand == -1 ]
  }
  return( SS_vec )
}
######

args <- commandArgs(trailingOnly = TRUE)
cn_file <- args[1] #~ 'test.txt'
exon_bed_file <- args[2] #~ '/home/sebastien/work_temp/id_card_temp/results/exons_lists/bed6/exon_up_siDNMT3b-siGL2_list.bed'
splicingSite <- args[3] #~ '3SS'

if ( cn_file == '-' ) {
  cn_tab <- read.table( file("stdin"), header=TRUE )
} else {
  cn_tab <- read.table( cn_file, header=TRUE )
}
exon_tab <- read.table( exon_bed_file, col.names=c( 'chr', 'start', 'end', 'exon_id', 'null', 'strand' ) )

merge_tab <- merge( cn_tab, exon_tab, by='exon_id', suffixes=c( '.cn', '' ) )


merge_tab[ 'SS' ] <- pos_SS( merge_tab, splicingSite )
merge_tab[ 'distToSS' ] <- ( merge_tab$pos - merge_tab$SS ) * merge_tab$strand

write.table( merge_tab[ , c( 'id', 'chr', 'pos', colnames( cn_tab )[ -c( 1,2,3,ncol( cn_tab ) ) ], 'distToSS', 'exon_id' ) ], file=stdout(), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE )
