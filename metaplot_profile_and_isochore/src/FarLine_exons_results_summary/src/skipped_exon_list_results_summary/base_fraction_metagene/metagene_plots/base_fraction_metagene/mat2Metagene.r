#!/usr/bin/Rscript

###########
mat2Metagene <- function ( mat, reduc_fc=function( xxx ) { mean( xxx, na.rm=TRUE ) } ) {
  # mat is the matrix to reduce in one metagene vector
  # reduc_fc is the function calculating one value per column of the input matrix
  apply( mat, 2, reduc_fc )
}
