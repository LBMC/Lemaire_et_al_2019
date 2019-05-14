#!/usr/bin/Rscript

## Script to build plot about one feature for constitutive exons (CE), alternatively spliced exons, and query list of exons
# Usage: Rscript exon_feature_plots.r <exons_features_table> <exon_list1,exon_list2,...> <feature1,feature2,...> --name <exon_name1,exon_name2,...>
# exons_features_table is the table of features for all the exons of FasterDB
# exon_list is a table with the columns 'id_gene', 'exon_pos', 'exons_flanquants'

# set environment
require( 'ggplot2' )
require( 'reshape2' )

all_args <- commandArgs(trailingOnly = F)
exon_feature_plots_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

## Some function
  # source( paste( exon_feature_plots_dir, 'val_retriever.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'order_conv.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'flanquing_exon_feat.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'intron_length.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'exon_tab_completer.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'feat_fig_builder.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'fix_strength_SS.r', sep='/' ) )
  source( paste( exon_feature_plots_dir, 'diff_strength_SS.r', sep='/' ) )


## Process arguments
args <- commandArgs(trailingOnly = TRUE)
feat_file <- args[1]
exon_list <- unlist( strsplit( args[2], split=',' ) )
out_dir <- args[3]
coding_exon_file <- args[4]

# ref_feat_tab_vec <- unlist( strsplit( args[5], split=',' ) )
# ref_name_vec <- unlist( strsplit( args[6], split=',' ) )
ref_feat_tab_vec <- c()
ref_name_vec <- c()
if ( '--ex-tab-ref' %in% args ) {
  index <- which( args == '--ex-tab-ref' )
  ref_feat_tab_vec <- unlist( strsplit( args[ index + 1 ], split=',' ) )
  ref_name_vec <- unlist( strsplit( args[ index + 2 ], split=',' ) )
}

ref_feat_tab_list <- list()
for ( xxx in ref_name_vec ) {
  tab_file <- ref_feat_tab_vec[ ref_name_vec == xxx ]
  ref_feat_tab_list[[ xxx ]] <- read.table( tab_file, header=TRUE, stringsAsFactors=FALSE )
}

# recover the name of the exons list
nom <- paste0( rep( "query", length( exon_list ) ), c( 1:length( exon_list ) ) )
if ( '--name' %in% args ) {
  nom <- unlist( strsplit( args[ which( args == '--name' ) + 1 ], split=',' ) )
}
nom_vec <- c( ref_name_vec, nom )

# recover the prefix of the output figures
fig_prefix <- 'test'
if ( '--fig-prefix' %in% args ) {
  fig_prefix <- args[ which( args == '--fig-prefix' ) + 1 ]
}

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}


## load the tables
# features table
feat_tab <- read.table( feat_file, header=TRUE, stringsAsFactors=FALSE, check.names=TRUE, na.strings=c( 'undef', 'exon_too_small' ) )
feat_tab[ , 'join_id' ] <- paste( feat_tab$GENE, feat_tab$EXON, sep='_' )

# exons table
exon_tab <- list()
for ( xxx in 1:length( exon_list ) ) {
  exon_tab[[ nom[ xxx ] ]] <- read.table( exon_list[ xxx ], header=TRUE, stringsAsFactors=FALSE )
  exon_tab[[ nom[ xxx ] ]][ , 'join_id' ] <- paste( exon_tab[[ nom[ xxx ] ]][ , 'id_gene' ], exon_tab[[ nom[ xxx ] ]][ , 'exon_pos' ], sep='_' )
  exon_tab[[ nom[ xxx ] ]] <- exon_tab_completer( exon_tab[[ nom[ xxx ] ]], feat_tab )
}

## compute the reference lists
final_data_list <- lapply( ref_feat_tab_list, function( xxx ) { feat_tab[ feat_tab$join_id %in% paste( xxx$GENE, xxx$EXON, sep='_' ), ] } )

# complete table of exon with the features (note: the list of exons is in fact a list of splicing events, exons with flanquing exons and specific introns )
final_data_list <- do.call( c, list( final_data_list, exon_tab ) )

## compute differential splicing site strength
final_data_list <- lapply( final_data_list, fix_strength_SS )
final_data_list <- lapply( final_data_list, diff_strength_SS )

# build the figure
feat_fig_builder( final_data_list, nom_vec, nom_query=c( nom ), out_dir, sign=fig_prefix, color_pallette=color_pallette )
