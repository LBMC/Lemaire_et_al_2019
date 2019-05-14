#!/usr/bin/Rscript

require( gplots, quietly=T, warn.conflicts=F )
require( RColorBrewer, quietly=T, warn.conflicts=F )

args <- commandArgs(trailingOnly = TRUE)

ratio_file_path <- args[ 1 ]
# ratio_file_path <- 'summary_mean_cov_ratio.tsv'
tab <- read.table( ratio_file_path, h=T, stringsAsFactors=F )
# print( str( tab ) )

if ( ncol( tab ) == 1 + 4 ) {
  tab <- cbind( tab, tab[ 5 ] )
}

out_png_pref <- args[ 2 ]
# out_png_pref <- 'bouh.png'

measure <- 'Mean'
flag <- '--measure'
if ( flag %in% args ) {
  measure <- args[ which( args == flag ) + 1 ]
}

mat_val <- as.matrix( tab[ tab$measure == measure, tail( colnames( tab ), n=-4 ) ] )
rownames( mat_val ) <- paste( tab[ tab$measure == measure, 'cond' ], tab[ tab$measure == measure, 'rep' ], tab[ tab$measure == measure, 'exp_des' ], sep='_' )
# str( mat_val )


#### heatmap example with color palette centered on 0
# heatmap.2(
#   mat_data,
#   cellnote = mat_data,  # same data set for cell labels
#   main = "Correlation", # heat map title
#   notecol="black",      # change font color of cell labels to black
#   density.info="none",  # turns off density plot inside color legend
#   trace="none",         # turns off trace lines inside the heat map
#   margins =c(12,9),     # widens margins around plot
#   col=my_palette,       # use on color palette defined earlier
#   breaks=col_breaks,    # enable color transition at specified limits
#   dendrogram="row",     # only draw a row dendrogram
#   Colv="NA"
# )            # turn off column clustering
#
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# col_breaks = c(seq(-1,0,length=100), # for red
# seq(0,0.8,length=100),  # for yellow
# seq(0.81,1,length=100)) # for green
####

sel_color_names <- as.character( floor( min( mat_val, na.rm=T ) * 100 ):ceiling( max( mat_val, na.rm=T ) * 100 ) )

nb_color <- 200
hm_col <- colorRampPalette( c( "blue", "white", "red" ) )( n = nb_color )
names( hm_col ) <- seq( -( nb_color / 2 ), ( nb_color / 2 ) - 1 )
hm_col <- hm_col[ head( sel_color_names, n=-1 ) ]

# col_breaks = c(
#   head( seq( min( c( mat_val, -0.1 ), na.rm=T ), 0, length=( nb_color / 2 ) + 1 ), n=-1 ),
#   0,
#   tail( seq( 0, max( c( mat_val, 0.1 ), na.rm=T ), length=( nb_color / 2 ) + 1 ), n=-1 )
# )

col_breaks = c(
  head( seq( -1, 0, length=( nb_color / 2 ) + 1 ), n=-1 ),
  0,
  tail( seq( 0, 1, length=( nb_color / 2 ) + 1 ), n=-1 )
)
names( col_breaks ) <- seq( -( nb_color / 2 ), nb_color / 2 )
col_breaks <- col_breaks[ sel_color_names ]


png( paste0( out_png_pref, '_', measure, '.png' ), heigh=1200, width=1200 )
# par( mar= c( 5, 5, 5, 5 ) )
heatmap.2(
  mat_val,
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins=c(5,30),     # widens margins around plot
  keysize=0.5,
  col=hm_col,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",
  Colv=NULL,
  Rowv=NULL,
  scale="none",
  na.col="#FFFFFF"
)
bouh <- dev.off()



####
