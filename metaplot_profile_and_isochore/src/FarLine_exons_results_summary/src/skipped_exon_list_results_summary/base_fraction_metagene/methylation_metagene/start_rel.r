#!/usr/bin/Rscript

all_args <- commandArgs(trailingOnly = F)
start_rel_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( start_rel_dir, 'reg_ext_window.r', sep='/' ) )


args <- commandArgs(trailingOnly = TRUE)
position <- as.numeric( args[1] )
window_size <- as.numeric( args[2] )
strand <- args[3]


tab_SS <- data.frame( start=c( position ), end=c( position ), strand=c( strand ), stringsAsFactors=FALSE )
tab_SS <- reg_ext_window( tab_SS, window_size )
write( tab_SS[ 1, 'start' ], file=stdout() )
