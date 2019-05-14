#!/usr/bin/Rscript

# require( gplots, quietly=T, warn.conflicts=F )
# require( RColorBrewer, quietly=T, warn.conflicts=F )

# require( ggplot2 )
# require( reshape2 )

args <- commandArgs(trailingOnly = TRUE)

summary_file_path <- args[ 1 ]
# ratio_file_path <- 'summary_mean_cov_ratio.tsv'
group_vec_arg <- args[ 2 ]
group_name_vec_arg <- args[ 3 ]
comp_list_arg <- args[ 4 ]

measure <- 'Mean' # 'Mean'
flag <- '--measure'
if ( flag %in% args ) {
  measure <- args[ which( args == flag ) + 1 ]
}


####
tab <- read.table( summary_file_path, h=T, stringsAsFactors=F )
tab[ 'sample' ] <- paste0( tab$cond, '_', tab$rep, '_', tab$exp_des )
# print( str( tab ) )

group_vec <- unlist( strsplit( group_vec_arg, split=',' ) )
group_name_vec <- unlist( strsplit( group_name_vec_arg, split=',' ) )
comp_list <- unlist( strsplit( comp_list_arg, split=',' ) )

write( paste( c( 'group', comp_list ), collapse='\t' ), file=stdout() )
for ( xxx in 1:length( group_vec ) ) {
  group <- read.table( group_vec[ xxx ], h=F, stringsAsFactors=F, sep='\t', fill=T )[ , 1 ]
  group <- group[ grep( pattern='^(#|$)', x=group, perl=T, invert=T ) ]
  group_name <- group_name_vec[ xxx ]
  sub_tab <- tab[ tab$exp_des %in% group , ]
  sub_tab <- sub_tab[ order( sub_tab$exp_des ), ]

  if ( nrow( sub_tab ) < 6 ) {
    write( paste0( '! Too few values for statistical test: ', group_name ), file=stderr() )
    next
  }

  pval_list <- list()
  out_line <- group_name
  for ( yyy in 1:length( comp_list ) ) {
    comp <- unlist( strsplit( comp_list[ yyy ], split='/' ) )
    set_list <- list()
    for ( zzz in comp ) {
      set_list[[ zzz ]] <- sub_tab[ grep( pattern=zzz, sub_tab$annot_list ), measure ]
    }
    pval_list[[ yyy ]] <- wilcox.test( set_list[[ 1 ]], set_list[[ 2 ]], paired=T )$p.value
    out_line <- paste0( out_line, '\t', pval_list[[ yyy ]] )
  }
  write( out_line, file=stdout() )
}



####
