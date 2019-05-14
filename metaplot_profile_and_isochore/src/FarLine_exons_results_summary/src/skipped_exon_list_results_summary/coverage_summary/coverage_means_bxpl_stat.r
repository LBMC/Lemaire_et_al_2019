#!/usr/bin/Rscript

require( reshape2 )
require( ggplot2 )

args <- commandArgs(trailingOnly = TRUE)
print( args )

mean_files <- unlist( strsplit( args[ 1 ], split=',' ) )
prefixes <- unlist( strsplit( args[ 2 ], split=',' ) )
out_dir <- args[ 3 ]
sign <- args[ 4 ]

color_pallette <- NULL
flag <- '--color-pallette'
if ( flag %in% args ) {
  index <- which( args == flag )
  print( index )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}
print( color_pallette )

mean_list <- list()
for ( xxx in 1:length( mean_files ) ) {
  mean_list[[ prefixes[ xxx ] ]] <- unlist( read.table( mean_files[ xxx ], h=F, stringsAsFactor=F ) )
  mean_list[[ prefixes[ xxx ] ]] <- as.numeric( mean_list[[ prefixes[ xxx ] ]] )
}

pval_file <- paste0( out_dir, '/pval_', sign, '.txt' ) #pval_wilcox.txt
write( 'statistical tests', file=pval_file )

# tests shapiro
write( 'shapiro tests', file=pval_file, append=T )
for ( xxx in 1:2 ) {
  if ( length( unique( mean_list[[ xxx ]] ) == 1 ) ) {
    pval <- 1
  } else {
    pval <- shapiro.test( mean_list[[ xxx ]] )[[ 'p.value' ]]
  }
  write( prefixes[ xxx ], file=pval_file, append=T )
  write( pval, file=pval_file, append=T )
  write( '', file=pval_file, append=T )
}

# comparison tests
write( paste0( prefixes[ 1 ], ' vs ', prefixes[ 2 ] ), file=pval_file, append=T )

# test Student/Welch
write( 'student test', file=pval_file, append=T )
stat_test <- t.test( mean_list[[ 1 ]], mean_list[[ 2 ]] )
pval <- stat_test[[ 'p.value' ]]
write( pval, file=pval_file, append=T )
write( '', file=pval_file, append=T )

# test wilcoxon
write( 'wilcoxon test', file=pval_file, append=T )
write( '-two.sided', file=pval_file, append=T )
stat_test <- wilcox.test( mean_list[[ 1 ]], mean_list[[ 2 ]], alternative='two.sided' )
pval <- stat_test[[ 'p.value' ]]
write( pval, file=pval_file, append=T )
write( '', file=pval_file, append=T )

write( '-greater', file=pval_file, append=T )
stat_test <- wilcox.test( mean_list[[ 1 ]], mean_list[[ 2 ]], alternative='greater' )
pval <- stat_test[[ 'p.value' ]]
write( pval, file=pval_file, append=T )
write( '', file=pval_file, append=T )

write( '-less', file=pval_file, append=T )
stat_test <- wilcox.test( mean_list[[ 1 ]], mean_list[[ 2 ]], alternative='less' )
pval <- stat_test[[ 'p.value' ]]
write( pval, file=pval_file, append=T )
write( '', file=pval_file, append=T )


## build plots
ggdata <- melt( mean_list )
ggdata$L1 <- factor( ggdata$L1, levels=prefixes )
fig <- ggplot( ggdata, aes( x=L1, y=value, fill=L1 ) )
fig <- fig + theme( axis.text=element_text( size=20 ), legend.position="none" )

if ( ! is.null( color_pallette ) ) {
  print( " in pallette" )
  fig <- fig + scale_fill_manual( values=color_pallette )
}

bxstat_list <- tapply( ggdata$value, ggdata$L1, boxplot.stats )
max_bxstat <- max( unlist( lapply( bxstat_list, function( xxx ) { xxx$stats[ 5 ] } ) ) )

fig_name <- paste0( out_dir, '/coverage_means_', sign, '_violin.png' )
fig_violin <- fig + geom_violin()
png( fig_name )
print(
  fig_violin
)
bouh <- dev.off()

fig_name <- paste0( out_dir, '/coverage_means_', sign, '_boxplot.png' )
fig_bxplt <- fig + geom_boxplot( outlier.shape=NA ) + coord_cartesian( ylim=c( 0, max_bxstat )*1.05 )
png( fig_name )
print(
  fig_bxplt
)
bouh <- dev.off()



####
