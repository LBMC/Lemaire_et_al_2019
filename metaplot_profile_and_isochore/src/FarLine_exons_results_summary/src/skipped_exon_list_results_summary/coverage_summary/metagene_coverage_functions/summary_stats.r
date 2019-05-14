#!/usr/bin/Rscript
f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  values <- as.numeric( unlist( strsplit( line, split='\t' ) ) )
  print( summary( values ) )
}
