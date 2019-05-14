#!/usr/bin/Rscript

# Load the libraries and custom functions
#~ library(rbamtools)
# set the environment
all_args <- commandArgs(trailingOnly = F)
coverage_means_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( coverage_means_dir, 'metagene_coverage_functions.r', sep='/' ) )


#########
## Recover the input arguments
harmonization <- FALSE
out_dir <- '.'
signature <- 'noS'
off_set <- 0

args <- commandArgs(trailingOnly = TRUE)
# print( args )
bw_file_arg <- args[ 1 ]
bed6_file_arg <- args[ 2 ]
out_file_arg <- args[ 3 ]

extra_args <- ''
if ( length( args > 3 ) ) {
  extra_args <- args[ 4:length( args ) ]
}


# load annotations
annot <- read.table( bed6_file_arg, header=FALSE, stringsAsFactors=FALSE )
annot <- annot[ ! grepl( '^#.*', annot[ , 1 ] ), ]
annot[ , 6 ] <- as.character( annot[ , 6 ] )
annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]
#########


script <- paste0( coverage_means_dir, '/metagene_coverage_functions/retrieve_bw_values_mean.py' )

annot_file <- tempfile()
write.table( annot, file=annot_file, quote=F, sep="\t", col.names=F, row.names=F)
#~ write.table( annot, file='./test.bed', quote=F, sep="\t", col.names=F, row.names=F)

out_file <- tempfile()

extra_args <- extra_args[ ! extra_args %in% c( '--stdout' ) ]
cmd=paste( 'python3', script, bw_file_arg, annot_file, out_file, paste( extra_args, collapse=' ' ), sep=' ' )
#~ if ( gene_concat ) {
#~   cmd=paste( cmd, '--group-concat', sep=' ' )
#~ }
# cmd=paste( script, bw_file, annot_file, out_file, nbins, sep=' ' )
message( cmd )
bouh <- system( cmd )

# recover out_file to write in stdout / out_file_arg if != '-'
con <- file( out_file, "r" )
lines <- paste( readLines( con, n=-1 ), sep='\n' )

out_file_real <- out_file_arg
message( out_file_real )
if ( out_file_arg == '-' & '--stdout' %in% args ) {
  out_file_real <- stdout()
}

write( lines, file=out_file_real )

####

