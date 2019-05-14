#!/usr/bin/Rscript

# Usage Rscript comp_deltaPSI_siRef_other.r <all_data_ref> <all_data_query> <tab_name1,tab_name2,...> [ --psi NUM ] [ --pval NUM ] [ --suf STR ] [ --out_dir PATH ] [ --corr | --antiC ] [ --no-filt-ref ] [ --no-filt-query ] [ --no-corr ] [ --no-enrich ] [ --no-counts ]


all_args <- commandArgs(trailingOnly = F)
comp_deltaPSI_siRef_other_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( comp_deltaPSI_siRef_other_dir, 'comp_deltaPSI_plot_function.r', sep='/' ) )
source( paste( comp_deltaPSI_siRef_other_dir, 'arg2value.r', sep='/' ) )


## Process args
args <- commandArgs(trailingOnly=TRUE)
all_data_ref <- args[1] #~ "/home/sebastien/FarLine_results_anaZip/FarLine_GSE52385_analysis"
all_data_query <- args[2]
tab_names <- unlist( strsplit( args[3], split=',' ) )

psi_thresh <- arg_value( '--psi', type_fc=as.numeric, default_value=0.1 )
pvalue_thresh <- arg_value( '--pval', type_fc=as.numeric, default_value=0.05 )
suffix <- arg_value( '--suf', default_value="noS" )

out_dir <- arg_value( '--out_dir', default_value="./")
dir.create( out_dir, showWarnings=FALSE )

tendency <- 'all' # allow to display all, correlated or anti-correlated exons
if ( '--corr' %in% args && '--antiC' %in% args ) {
  stop( '--corr and --antiC options cannot exist in the same time, write nothing to take all exons in account.' )
} else if ( '--corr' %in% args ) {
  tendency <- 'corr'
} else if ( '--antiC' %in% args ) {
  tendency <- 'antiC'
}

common_only <- FALSE # display only common filtered exons
if ( '--common-only' %in% args ) {
  common_only <- TRUE
}

filt_query <- TRUE # filter for significative exon in query analysis
if ( '--no-filt-query' %in% args ) {
  filt_query <- FALSE
}
filt_ref <- TRUE # filter for significative exon in ref analysis
if ( '--no-filt-ref' %in% args ) {
  filt_ref <- FALSE
}

pval_highl <- FALSE # display counts for common sig exons ( useful if --common-only and --no-filt* is used )
if ( '--pval-highlight' %in% args ) {
  pval_highl <- TRUE
}

correlation <- TRUE # display test results for correlation
if ( '--no-corr' %in% args ) {
  correlation <- FALSE
}

enrichment <- TRUE # display test results for enrichments
if ( '--no-enrich' %in% args ) {
  enrichment <- FALSE
}

comptage <- TRUE # display counts of exons ( specific exons, common exons )
if ( '--no-counts' %in% args ) {
  comptage <- FALSE
}

## Process data
siRef <- read.table( all_data_ref, header=TRUE, stringsAsFactors=FALSE, row.names=c( "id_event" ) )
cust <- read.table( all_data_query, header=TRUE, stringsAsFactors=FALSE, sep='\t', row.names=c( "id_event" ) )

## build the plot
fig_name <- paste( out_dir, "/", suffix, ".png", sep='' )
# message( fig_name )
png( fig_name, width=1200, height=1200, pointsize=24 )
comp_deltaPSI_plot( siRef, cust, first_name=tab_names[ 1 ], second_name=tab_names[ 2 ], psi_thresh=psi_thresh, pvalue_thresh=pvalue_thresh, tendency=tendency, common_only=common_only, filt_first=filt_ref, filt_second=filt_query, pval_highl=pval_highl, correlation=correlation, enrichment=enrichment, comptage=comptage )
bouh <- dev.off()

write( fig_name, file=stdout() )
