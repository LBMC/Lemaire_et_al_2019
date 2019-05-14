#!/usr/bin/Rscript

rnafold2MFEtable <- function( filepath ) {
  mfe_df <- list()

  con <- file(filepath, "r")

  ## parse the file of RNAFold output
  nnn <- 0
  while ( TRUE ) {
    line <- readLines(con, n = 1)

    if ( length(line) == 0 ) {
      break
    }

    nnn <- nnn + 1
    if ( startsWith( x=line, prefix='>' ) ) {
      nnn <- 1
      nom <- paste( tail( strsplit( line, split='' )[[ 1 ]], n=-1 ), collapse='' )
    }

    if ( nnn == 3 ) {
      mfe_df[[ nom ]] <- head( strsplit( tail( strsplit( line, split='(', fixed=TRUE )[[ 1 ]], n=1 ), split=')', fixed=TRUE )[[ 1 ]], n=1 )
    }

    # print( line )
  }

  close(con)

  ## restructure the recovered data
  mfe_df <- melt( mfe_df )
  colnames( mfe_df ) <- c( 'MFE', 'id' )
  mfe_df[ , 'MFE' ] <- as.numeric( as.character( mfe_df[ , 'MFE' ] ) )

  return( mfe_df )
}
