#!/usr/bin/Rscript

#### compute summary over the table of coverages
summary_table <- function( cov_df, repl ) {
    tmp <- lapply( cov_df, summary )
    mm <- melt( do.call( rbind, tmp ) )
    colnames( mm )[ colnames( mm ) == "X2" ] <- "measure"
    mm <- data.frame( mm, replicate = rep( repl, dim( mm )[1] ) )
    return( mm )
}
