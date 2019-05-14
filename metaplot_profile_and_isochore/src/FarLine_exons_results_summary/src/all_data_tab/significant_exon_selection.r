#!/usr/bin/Rscript

##############
## Written by Sébastien Lemaire
## on 06/07/2014
##############

## Script to select the significative exons from a table under FarLine format ( default is p-value <= 0.05, deltaPSI >= 0.1 )

## Usage: [Rscript] significant_exon_selection.r input_tab [ --pval NUM, --delta NUM ]
# --pval is to set the maximum p-value
# --delta is to set the minimum deltaPSI

# the columns "valeur.abs.deltaPSI" and "pvalue_corrigee_glm", or"pvalue_corrigee_Fisher", have to be present
# !!! if columns named "pvalue_corrigee_glm" is present, this column is considered instead of "pvalue_corrigee_Fisher" if also present
# the output table is written in stdout

## Start of the script
# treat the arguments
args <- commandArgs( trailingOnly=TRUE )

input_table <- args[ 1 ]


thres_pval <- 0.05
if ( '--pval' %in% args ) {
	thres_pval <- as.numeric( args[ which( args == '--pval' ) + 1 ] )
}

thres_delta <- 0.1
if ( '--delta' %in% args ) {
	thres_delta <- as.numeric( args[ which( args == '--delta' ) + 1 ] )
}

# load the input table
if ( input_table == '-') {
  # message( 'Input: stdin' )
  tab <- read.table(file("stdin"), sep = '\t', header = TRUE, stringsAsFactors=FALSE )
} else {
  #~ message( paste( 'Input: ', args[1], sep='') )
  tab <- read.table(input_table, sep = '\t', header = TRUE, stringsAsFactors=FALSE )
}

#colonnes <- c( "gene_symbol", "exon_skipped", "coordonnees", "exons_flanquants", "deltaPSI", "valeur.abs.deltaPSI",)

# determine presence of "pvalue_corrigee_glm" column
if ( "pvalue_corrigee_glm" %in% colnames( tab ) ) {
	pcol <- "pvalue_corrigee_glm"
} else {
	pcol <- "pvalue_corrigee_Fisher"
}

# select the significative exons
tab_sig <- tab[ tab[ , pcol ] <= thres_pval & !is.na( tab[ , pcol ] ) & tab[ , "valeur.abs.deltaPSI" ] >= thres_delta & !is.na( tab[ , "valeur.abs.deltaPSI" ] ), ]

# print to stdout the results
write.table( x=tab_sig, file=stdout(), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t' )
