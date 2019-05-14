#!/usr/bin/Rscript

# require( gplots, quietly=T, warn.conflicts=F )
# require( RColorBrewer, quietly=T, warn.conflicts=F )

require( ggplot2 )
require( reshape2 )

args <- commandArgs(trailingOnly = TRUE)

ratio_file_path <- args[ 1 ]
# ratio_file_path <- 'summary_mean_cov_ratio.tsv'
group_vec_arg <- args[ 2 ]
group_name_vec_arg <- args[ 3 ]
out_png_pref <- args[ 4 ]
# out_png <- 'bouh.png'

measure <- 'Mean' # 'Mean'
flag <- '--measure'
if ( flag %in% args ) {
  measure <- args[ which( args == flag ) + 1 ]
}


####
tab <- read.table( ratio_file_path, h=T, stringsAsFactors=F )
# print( str( tab ) )

group_vec <- unlist( strsplit( group_vec_arg, split=',' ) )
group_name_vec <- unlist( strsplit( group_name_vec_arg, split=',' ) )

ratio_cols <- 5:ncol( tab )
for ( xxx in ratio_cols ) {
  tab[ , xxx ] <- as.numeric( tab[ , xxx ] )
}
tab[ 'group' ] <- NA

group_list <- list()
for ( xxx in 1:length( group_vec ) ) {
  gname <- group_name_vec[ xxx ]
  group_list[[ gname ]] <- unlist( read.table( group_vec[ xxx ], h=F, stringsAsFactors=F ) )

  tab[ tab[ , 'exp_des' ] %in% group_list[[ gname ]], 'group' ] <- gname
}

tab <- tab[ ! is.na( tab$group ), ]
tab[ , 'group' ] <- factor( tab[ , 'group' ], levels=group_name_vec )
# print( str( tab ) )

####

# options( width=250 )
#
# stab <- tab[ tab$measure == measure & tab$group == "H3K79me2", ]
#
# for ( g_name in group_name_vec ) {
#   stab <- tab[ tab$measure == measure & tab$group == g_name, ]
#   if ( nrow( stab ) >= 3 ){
#     print( g_name )
#     print( shapiro.test( stab[ , 5 ] )[[ 'p.value' ]] )
#     print( '' )
#   }
# }

####

for ( ratio_col_plot in ratio_cols ) {
  # ratio_col_plot <- 7
  gtab <- tab[ tab$measure == measure ,c( 1:4, ratio_col_plot, ncol( tab ) ) ]
  colnames( gtab ) <- c( colnames( tab )[ 1:4 ], 'ratio', 'group' )
  fig <- ggplot( gtab, aes( x=group, y=ratio ) ) + geom_boxplot() + ylim( -1, 1 )
  # str( gtab )

  png_name <- paste0( out_png_pref, '_', measure, '_', colnames( tab )[ ratio_col_plot ], '.png' )
  # png_name <- paste0( out_png_pref, '_', measure, '_', as.character( ratio_col_plot ), '.png' )
  message( png_name )
  png( png_name )
  print( fig + theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) )
  bouh <- dev.off()

  # stop()
}





####
