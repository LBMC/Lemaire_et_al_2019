#!/usr/bin/Rscript

# Load the libraries and custom functions
#~ library(rbamtools)
require(bigWig, quietly=TRUE, warn.conflicts=FALSE )
require(ggplot2, quietly=TRUE, warn.conflicts=FALSE )
require(reshape, quietly=TRUE, warn.conflicts=FALSE )
require( purrr, quietly=TRUE, warn.conflicts=FALSE )
# require(gplots, quietly=TRUE, warn.conflicts=FALSE )

# set the environment
all_args <- commandArgs(trailingOnly = F)
metagene_coverage_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( metagene_coverage_dir, 'metagene_coverage_functions.r', sep='/' ) )


#########
# Recover the input arguments
harmonization <- FALSE
out_dir <- '.'
signature <- 'noS'

args <- commandArgs(trailingOnly = TRUE)
# print( args )
bw_files <- unlist( strsplit( args[ which( args == '-bw' ) + 1 ], split=',' ) ) #~ list of the bw files (~ sample) to consider
conditions <- unlist( strsplit( args[ which( args == '-cond' ) + 1 ], split=',' ) ) #~ indicate for each bw_files
replicates <- unlist( strsplit( args[ which( args == '-rep' ) + 1 ], split=',' ) ) #~ replicate for each bw_file

annot_file_vec <- unlist( strsplit( args[ which( args == '-annot' ) + 1 ], split=',' ) ) #~ bed file of the annotations to consider
prefixes <- unlist( strsplit( args[ which( args == '-prefixes' ) + 1 ], split=',' ) )
comp_pair <- unlist( strsplit( args[ which( args == '-comp_pair' ) + 1 ], split=',' ) )
message( paste( 'comp_pair', paste( comp_pair, collapse=' ' ) ) )



# optional arguments
if ( '-sign' %in% args ) {
index <- which( args == '-sign' )
signature <- args[ index + 1 ]
}
if ( '-off_set' %in% args ) {
index <- which( args == '-off_set' )
off_set <- as.numeric( args[ index + 1 ] )
}
if ( '-out_dir' %in% args ) {
index <- which( args == '-out_dir' )
out_dir <- args[ index + 1 ]
}
ylims_list=list( 'mean'=c( NA, NA ), 'median'=c( NA, NA ) )
if ( '--ylims_mean' %in% args ) {
index <- which( args == '--ylims_mean' )
ylims_list[[ 'mean' ]] <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

ylims_median=c( NA, NA )
if ( '--ylims_median' %in% args ) {
index <- which( args == '--ylims_median' )
ylims_list[[ 'median' ]] <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

if ( '--scaleHarm' %in% args ) {
  harmonization <- TRUE
  ylims_list <- list( mean=c( -1, 1 ), median=c( -1, 1 ) )
}


# create the output folders of the images
list_out_dir <- out_dirs_creator( out_dir, all=FALSE )
#########


#########
## set the environment
# load the sequencing depth information
df_expdesign <- data.frame( file=bw_files, condition=conditions, replicate=replicates, row.names=seq( 1, length( bw_files ) ) )
df_expdesign$file <- as.character( df_expdesign$file )
vec_seqDepth <- sapply( df_expdesign[[ 'file' ]], bw_seqDepth_retriever )
seqDepth_df <- cbind( df_expdesign, sequencing_depth=vec_seqDepth )
seqDepth_df[ 'normFactor' ] <- seqDepth_df[ ,'sequencing_depth' ]/min( seqDepth_df[ ,'sequencing_depth' ] )
message( '> Experiment design')
print( seqDepth_df )

# write the data frame of experiment design
write.table( seqDepth_df, file=paste( out_dir, 'seqDepth.tsv', sep='/' ), sep='\t', quote=FALSE, row.names=FALSE )


# load annotations
names( annot_file_vec ) <- prefixes
annot_list <- lapply( annot_file_vec, function( annot_file ) {
  annot <- read.table( annot_file, header=FALSE, stringsAsFactors=FALSE )
  annot[ , 6 ] <- as.character( annot[ , 6 ] )
  annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]
  return( annot )
} )
#########


#########
## retrive the coverages and compute the resuming values
fc_vec <- c( 'mean', 'median' )
list_total_mm <- list()
for ( prefix in prefixes ) {
  # message( prefix )
  list_mm <- list()

  for ( repl in levels( as.factor( seqDepth_df[ seqDepth_df$condition %in% comp_pair,'replicate' ] ) ) ) {
    # message( repl )
    bw_file_list <- list()
    normFactor_list <- list()

    for ( cond in comp_pair ) {
      # message( cond )
      bw_file_list[[ cond ]] <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]
      normFactor_list[[ cond ]] <- seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ]
    }
      # retrieve coverage from bam file for all annotations and apply normalization
      diff_cov_df <- stranded_diffCov( bw_file_list, annot_list[[ prefix ]], off_set=off_set, normFactor_list, comp_pair, harmonization=harmonization )

    # resume the coverages per position
    for ( type_fc in fc_vec ) {
      list_mm[[ type_fc ]][[ repl ]] <- apply( diff_cov_df, 2, eval( parse( text=type_fc ) ) )
    }
  }
  # message( 'sub' )
  # print( str( list_mm[[ 'mean' ]] ) )

  # resume the replicates the median applied on them for each position
  for ( type_fc in fc_vec ) {
    list_total_mm[[ type_fc ]][[ paste( rev( comp_pair ), collapse='-' ) ]][[ prefix ]] <- apply( do.call( rbind, list_mm[[ type_fc ]] ), 2, median )
  }
}

# structure the resuming value per annotation list in one table
list_total_mm <- list_total_mm %>% modify_depth( .depth=2, .f=function ( xxx ) { do.call( rbind, xxx ) } )

# message( 'total' )
# print( str( list_total_mm ) )
# saveRDS( list_total_mm, 'bouh.RDS' )
#########


#########
## build the figures
for ( type_fc in c( 'mean', 'median' ) ) {
    # curve
  total_mm <- list_total_mm[[ type_fc ]] %>% melt()
  colnames( total_mm ) <- c( 'annot', 'X1', 'value', 'condition')
  total_mm[ , 'annot' ] <- factor( total_mm[ , 'annot' ], levels=prefixes )
  # message( 'melt' )
  # print( str( total_mm ) )
  metagene_measure_builder( paste( list_out_dir[[ paste( 'metagene_', type_fc, sep='' ) ]], '/metagene_', type_fc, '_', signature, '.png' , sep='' ), total_mm, ylims=ylims_list[[ type_fc ]], colour=total_mm$annot, color_pallette=color_pallette )

    # heatmap
  for ( cond in names( list_total_mm[[ type_fc ]] ) ) {
    heatmap_builder( paste( list_out_dir[[ 'heatmapSimple' ]], '/heatmapSimple_', type_fc, '_', cond, '_', signature, '.png', sep='' ), list_total_mm[[ type_fc ]][[ cond ]] )
  }
}
# stop( 'bouh' )
#########
