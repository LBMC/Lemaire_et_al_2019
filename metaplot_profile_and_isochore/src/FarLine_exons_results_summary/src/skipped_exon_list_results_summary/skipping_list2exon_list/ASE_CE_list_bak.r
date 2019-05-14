#!/usr/bin/Rscript

## Process argumentsF
args <- commandArgs(trailingOnly = TRUE)
feat_file <- args[1]
coding_exon_file <- args[2]
exon_type <- args[3]

# load tables
feat_tab <- read.table( feat_file,
  header=TRUE,
  stringsAsFactors=FALSE,
  check.names=TRUE,
  na.strings=c( 'undef', 'exon_too_small' )
  )
feat_tab[ , 'join_id' ] <- paste( feat_tab[ , 'GENE' ], feat_tab[ , 'EXON' ], sep='_' )

coding_exons <- read.table( coding_exon_file,
  header=TRUE,
  stringsAsFactors=FALSE,
  sep=';',
  quote='"'
  )
coding_exons[ 'join_id' ] <- paste( coding_exons$id_gene, coding_exons$pos_sur_gene, sep='_' )

# retrieve the ids
if ( exon_type == 'ASE' ) {
  exon_ids <- coding_exons$join_id[ coding_exons$exon_types == 'ACE' ]
} else if ( exon_type == 'CE' ) {
  exon_ids <- coding_exons$join_id[ coding_exons$exon_types == 'CCE' ]
}

# filter feat_tab for the ids
exontab <- feat_tab[ feat_tab$join_id %in% exon_ids,  ]

# write the tables
write.table( x=exontab, file=stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
