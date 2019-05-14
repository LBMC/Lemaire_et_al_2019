#!/usr/bin/Rscript

order_conv_dir <- dirname(sys.frame(1)$ofile)
source( paste( order_conv_dir, 'val_retriever.r', sep='/' ) )

# give the order allowing to reorder a vector according to another one containing the same elements with the same occurences
order_conv <- function ( query, ref ) {
  order( query )[ order( order( ref ) ) ]
}
