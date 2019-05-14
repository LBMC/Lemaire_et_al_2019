#!/usr/bin/Rscript

## Script to gather all the data from exon_skipping FarLine analysis in one table
# Usage: [Rscript] all_data_computer.r <analyse_stat directory> <output_directory>

# the "analyse_stat" directory is a directory given by the statistical analysis of FarLine on the RNA-Seq data.
# the output file "all_data.tab" is written in the output directory, you have to rename it by yourself then.

## Main Script
message( 'Start of script.' )

# Recover command line arguments
args <- commandArgs(trailingOnly = TRUE)
comparaison_dir <- args[ 1 ]
out_dir <- args[ 2 ]

# general variables
skipping_folder <- paste( comparaison_dir, 'exon_skipping/tmp', sep='/' )

all_count_tab <- read.table( paste( skipping_folder, "all_count.tab", sep='/' ), header=TRUE, stringsAsFactors=FALSE, sep='\t' )
tmp_event_tab <- read.table( paste( skipping_folder, "tmp_event.tab", sep='/' ), header=TRUE, stringsAsFactors=FALSE, sep='\t' )
tmp_final_data_tab <- read.table( paste( skipping_folder, "tmp_final_data.tab", sep='/' ), header=TRUE, stringsAsFactors=FALSE, sep='\t' )

# order the tables for "id_event" field
all_count_tab <- all_count_tab[ order( all_count_tab$id_event ), ]
tmp_event_tab <- tmp_event_tab[ order( tmp_event_tab$id_event ), ]
tmp_final_data_tab <- tmp_final_data_tab[ order( tmp_final_data_tab$id_event ), ]

#~ print( str( all_count_tab ) )
#~ print( str( tmp_event_tab ) )
#~ print( str( tmp_final_data_tab ) )
#~ print( all( all_count_tab$id_event == tmp_event_tab$id_event ) )
#~ print( all( all_count_tab$id_event == tmp_final_data_tab$id_event ) )

# compute the PSI values from the counts
psi_tab <- all_count_tab[ , grep( 'inclu', colnames( all_count_tab ) ) ] / ( all_count_tab[ , grep( 'inclu', colnames( all_count_tab ) ) ] + all_count_tab[ , grep( 'exclu', colnames( all_count_tab ) ) ] )
colnames( psi_tab ) <- gsub( 'inclu', 'psi', colnames( all_count_tab )[ grep( 'inclu', colnames( all_count_tab ) ) ] )

# gather in one table
complete_tab <- cbind( tmp_final_data_tab, id_gene=tmp_event_tab$id_gene, psi_tab, all_count_tab[ , ! colnames( all_count_tab ) %in% c( "id_event" ) ] )

#~ print( str( complete_tab ) )

# write the complete_tab in a file
# write.table( complete_tab, file=paste( out_dir, 'all_data.tab', sep='/' ), sep='\t', quote=FALSE, row.names=FALSE )
write.table( complete_tab, file=stdout(), sep='\t', quote=FALSE, row.names=FALSE )

message( 'End of script.' )
