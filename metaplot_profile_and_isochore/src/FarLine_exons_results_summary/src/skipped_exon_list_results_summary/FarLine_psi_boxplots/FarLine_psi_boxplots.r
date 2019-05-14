#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

all_args <- commandArgs(trailingOnly = F)
FarLine_psi_boxplots_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

## Some function
source( paste( FarLine_psi_boxplots_dir, 'exonTab_FarLineTab_joiner.r', sep='/' ) )
source( paste( FarLine_psi_boxplots_dir, 'FasterDBTab_FarLineTab_joiner.r', sep='/' ) )
source( paste( FarLine_psi_boxplots_dir, 'FarLine_psi_summary.r', sep='/' ) )
source( paste( FarLine_psi_boxplots_dir, 'FarLine_psi_plots.r', sep='/' ) )
source( paste( FarLine_psi_boxplots_dir, 'FarLine_psi_reparts.r', sep='/' ) )

## process arguments and load tables
args <- commandArgs(trailingOnly = TRUE)
exons_table_vec <- unlist( strsplit( args[1], split=',' ) )
FarLine_all_data_tab <- read.table( args[2], header=TRUE, stringsAsFactors=FALSE )
feat_tab <- read.table( args[3], header=TRUE, stringsAsFactors=FALSE )
coding_exon_tab <- read.table( args[4], header=TRUE, stringsAsFactors=FALSE, sep=';', quote='"' )
out_dir <- args[5]

# ref_feature_tab <- unlist( strsplit( args[6], split=',' ) )
# ref_name <- unlist( strsplit( args[7], split=',' ) )
ref_feat_tab <- c()
ref_name <- c()
if ( '--ex-tab-ref' %in% args ) {
  index <- which( args == '--ex-tab-ref' )
  ref_feat_tab <- unlist( strsplit( args[ index + 1 ], split=',' ) )
  ref_name <- unlist( strsplit( args[ index + 2 ], split=',' ) )
}

ref_feat_tab_list <- list()
for ( xxx in ref_name ) {
  tab_file <- ref_feature_tab[ ref_name == xxx ]
  ref_feat_tab_list[[ xxx ]] <- read.table( tab_file, header=TRUE, stringsAsFactors=FALSE )
}

# recover the name of the exons list
nom <- c( "query" )
if ( '--name' %in% args ) {
  nom <- unlist( strsplit( args[ which( args == '--name' ) + 1 ], split=',' ) )
}
nom_vec <- c( ref_name, nom )

# recover the prefix of the output figure
fig_prefix <- "test"
if ( '--fig-prefix' %in% args ) {
  fig_prefix <- args[ which( args == '--fig-prefix' ) + 1 ]
}

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}


## merge/join the tables
# extraction FarLine data for ASE and CE separately
feat_complete_tab <- FasterDBTab_FarLineTab_joiner( feat_tab, FarLine_all_data_tab )
list_complete_tab <- lapply( ref_feat_tab_list, function( xxx ) { feat_complete_tab[ feat_complete_tab$join_id %in% paste( xxx$GENE, xxx$EXON ), ] } )

# add the query exon table FarLine table
for ( xxx in c( 1:length( nom ) ) ) {
  exons_table <- read.table( exons_table_vec[ xxx ], header=TRUE, stringsAsFactors=FALSE )
  list_complete_tab[[ nom[ xxx ] ]] <- exonTab_FarLineTab_joiner( exons_table, FarLine_all_data_tab )
}

## compute the median of the psi per condition
list_psi_median_tab <- lapply( list_complete_tab, FarLine_psi_summary, reduc_fc=median, na.rm=TRUE )

## compute the output suffix
suffix <- paste( colnames( list_psi_median_tab[[ 1 ]] ), collapse='-' )

# build the boxplots, violin plots
draw_quantiles=c( 0.25, 0.5, 0.75 )
FarLine_psi_plots( list_psi_median_tab, nom_vec, out_dir, sign=fig_prefix, suffix=suffix, draw_quantiles=draw_quantiles, color_pallette=color_pallette )

## compute the repartitions
list_psi_repart_tab <- lapply( list_psi_median_tab,
  function( xxx ) {
    new_xxx <- lapply( xxx, function( yyy ) {
      new_yyy <- data.frame( psi=yyy, repart=order( order( yyy ) )/length( yyy ) )
      return( new_yyy )
    } )
    return( new_xxx )
  } )

## build the repartition plots
FarLine_psi_reparts( list_psi_repart_tab, nom_vec, out_dir, sign=fig_prefix, color_pallette=color_pallette )
