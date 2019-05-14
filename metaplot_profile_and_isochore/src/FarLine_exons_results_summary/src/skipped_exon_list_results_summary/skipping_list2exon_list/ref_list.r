#!/usr/bin/Rscript

## Script to recover sub part of the FasterDB features table from a list of exon identifiers (format: geneId_exonIdinGene, with headers)

options(scipen=999)

## Process argumentsF
args <- commandArgs(trailingOnly = TRUE)
feat_file <- args[1]
ref_exons_list <- args[2]
# coding_exon_file <- args[2]
# exon_type <- args[3]

# load tables
feat_tab <- read.table( feat_file,
  header=TRUE,
  stringsAsFactors=FALSE,
  check.names=TRUE,
  na.strings=c( 'undef', 'exon_too_small' )
  )
feat_tab_join_id <- paste( feat_tab[ , 'GENE' ], feat_tab[ , 'EXON' ], sep='_' )

ref_exons <- read.table( ref_exons_list,
  header=TRUE,
  stringsAsFactors=FALSE,
  )[ , 1 ]

# retrieve the ids
# if ( exon_type == 'ASE' ) {
#   exon_ids <- ref_exons$join_id[ ref_exons$exon_types == 'ACE' ]
# } else if ( exon_type == 'CE' ) {
#   exon_ids <- ref_exons$join_id[ ref_exons$exon_types == 'CCE' ]
# }

# filter feat_tab for the ids
exontab <- feat_tab[ feat_tab_join_id %in% ref_exons,  ]
# exontab <- feat_tab[ feat_tab$join_id %in% exon_ids,  ]

# write the tables
write.table( x=exontab, file=stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
