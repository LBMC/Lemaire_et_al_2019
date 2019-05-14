#!/usr/bin/Rscript

## set environment
options(scipen=999)

# find directory of the script
all_args <- commandArgs(trailingOnly = F)
reg_ss_coord_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

# remove scrientific notation
options(scipen=999)

## Some function
# source( paste( reg_ss_coord_dir, 'ss_reg.r', sep='/' ) )
source( paste( reg_ss_coord_dir, 'reg_ext_window.r', sep='/' ) )


## process arguments
args <- commandArgs(trailingOnly = TRUE)

tab_file <- args[1] #~ './results/exons_lists/ASE_FDB_list.tsv'
window_size <- as.numeric( args[2] )


## load the table
tab <- read.table( tab_file, header=TRUE, stringsAsFactors=FALSE )

## process the table
coord_table <- do.call( rbind, strsplit( tab$coord, split=':' ) )
border_table <- do.call( rbind, strsplit( coord_table[ , 2 ], split='-' ) )
coord_table <- as.data.frame( cbind( coord_table[ , 1 ], border_table, tab$strand ), stringsAsFactors=FALSE )
coord_table[ , 2 ] <- as.numeric( coord_table[ , 2 ] )
coord_table[ , 3 ] <- as.numeric( coord_table[ , 3 ] )
colnames( coord_table ) <- c( 'ref', 'start', 'end', 'strand' )

## extend the splicing site regions for window centered on borders
coord_table <- reg_ext_window( coord_table, window_size )


## write table of splicing sites
  tab_ex <- tab
  tab_ex$coord <- paste( coord_table$start, coord_table$end, sep='-' )
  tab_ex$coord <- paste( coord_table$ref, tab_ex$coord, sep=':' )
  write.table( x=tab_ex, file=stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
