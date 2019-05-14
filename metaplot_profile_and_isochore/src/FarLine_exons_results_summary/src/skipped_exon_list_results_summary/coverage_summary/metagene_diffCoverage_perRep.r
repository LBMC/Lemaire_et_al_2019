#!/usr/bin/Rscript

# Load the libraries and custom functions
#~ library(rbamtools)
require(bigWig, quietly=TRUE, warn.conflicts=FALSE )
require(ggplot2, quietly=TRUE, warn.conflicts=FALSE )
require(reshape, quietly=TRUE, warn.conflicts=FALSE )
require(gplots, quietly=TRUE, warn.conflicts=FALSE )

# set the environment
all_args <- commandArgs(trailingOnly = F)
metagene_coverage_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
  source( paste( metagene_coverage_dir, 'metagene_coverage_functions.r', sep='/' ) )

# Recover the input arguments
harmonization <- FALSE
out_dir <- '.'
signature <- 'noS'

args <- commandArgs(trailingOnly = TRUE)
bw_files <- unlist( strsplit( args[ which( args == '-bw' ) + 1 ], split=',' ) ) #~ list of the bw files (~ sample) to consider
conditions <- unlist( strsplit( args[ which( args == '-cond' ) + 1 ], split=',' ) ) #~ indicate for each bw_files
replicates <- unlist( strsplit( args[ which( args == '-rep' ) + 1 ], split=',' ) ) #~ replicate for each bw_file

annot_file <- args[ which( args == '-annot' ) + 1] #~ bed file of the annotations to consider
comp_pair <- unlist( strsplit( args[ which( args == '-comp_pair' ) + 1 ], split=',' ) )


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
if ( '-scaleHarm' %in% args ) {
harmonization <- TRUE
}

ylims_mean=c( NULL, NULL )
if ( '-ylims_mean' %in% args ) {
index <- which( args == '-ylims_mean' )
ylims_mean <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

ylims_median=c( NULL, NULL )
if ( '-ylims_median' %in% args ) {
index <- which( args == '-ylims_median' )
ylims_median <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

# create the output folders of the images
list_out_dir <- out_dirs_creator( out_dir )


#### Create the figures
# load the sequencing depth information
df_expdesign <- data.frame( file=bw_files, condition=conditions, replicate=replicates, row.names=seq( 1, length( bw_files ) ) )
df_expdesign$file <- as.character( df_expdesign$file )
vec_seqDepth <- sapply( df_expdesign[[ 'file' ]], bw_seqDepth_retriever )
seqDepth_df <- cbind( df_expdesign, sequencing_depth=vec_seqDepth )
seqDepth_df[ 'normFactor' ] <- seqDepth_df[ ,'sequencing_depth' ]/min( seqDepth_df[ ,'sequencing_depth' ] )
message( '> Experiment design')
print( seqDepth_df )


# load annotations
annot <- read.table( annot_file, header=FALSE, stringsAsFactors=FALSE )
annot[ , 6 ] <- as.character( annot[ , 6 ] )
# annot <- cbind( annot, strand=do.call( rbind, strsplit( annot[ ,5 ], split='ppp' ) )[ ,2 ] )
annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]


list_total_mm <- list()

  list_mm <- list()
  for ( xxx in levels( as.factor( seqDepth_df[ seqDepth_df$condition %in% comp_pair,'replicate' ] ) ) ) {
    bw_file_list <- list()
    cov_df_list <- list()
    normFactor_list <- list()
    for (cond in comp_pair) {
      # bw_pattern <- paste( cond, 'n', xxx, ".*_treat_pileup.bw$", sep='' )
      # print( bw_pattern )
      # bw_name <- list.files( path = bw_dir, pattern = bw_pattern )
      # bw_file_list[[ cond ]] <- paste( bw_dir, bw_name, sep='/' )
      bw_file_list[[ cond ]] <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == xxx ]

      normFactor_list[[ cond ]] <- seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==xxx, 'normFactor' ]

    }
    # print( 'BigWig in process' )
    # print( bw_file_list )

    # retrieve coverage from bam file for all annotations and apply normalization
    diff_cov_df <- stranded_diffCov( bw_file_list, annot, off_set=off_set, normFactor_list, comp_pair, harmonization=harmonization )

    # apply normalizatio    # per-base summary of coverage, one curve per statistics
    list_mm[[ xxx ]] <- summary_table( diff_cov_df, xxx )

    # Build heatmap
    heatmap_builder( paste( list_out_dir[[ 'heatmapSimple' ]], '/heatmapSimple_diffCov_', comp_pair[ 1 ], 'n', xxx, '_', comp_pair[ 2 ], 'n', xxx, '_', signature, '.png', sep='' ), diff_cov_df )

    #~ png( filename = paste( './', heatmapNormalized_dir, '/heatmapNormalized_', cell_type, cond, 'n', xxx, '_', signature, '.png', sep='' ), width = 1920, height = 1200, pointsize = 48 )
    #~ temp_sums <- rowSums( cov_df )
    #~ temp_sums <- temp_sums
    #~ temp_sums[ temp_sums == 0 ] <- 1
    #~ heatmap( sweep( as.matrix( cov_df ), 1, temp_sums, '/' ), cexCol=0.2, Colv=NA )
    #~ dev.off()

  }

  # plot the summaries
  # list_total_mm[[ cond ]] <- do.call( rbind, list_mm )
  list_total_mm <- do.call( rbind, list_mm )

  metagene_summary_builder( paste( list_out_dir[[ 'metagene_summary' ]], '/metagene_summary_diffCov_', comp_pair[ 1 ], '-', comp_pair[ 2 ], '_', signature, '.png' , sep='' ), list_total_mm )

# plot the means
# total_mm_mean <- measure_extractor( list_total_mm, "Mean", diff=TRUE )
total_mm_mean <- measure_extractor( list_mm, "Mean", diff=TRUE )

metagene_measure_builder( paste( list_out_dir[[ 'metagene_mean' ]], '/metagene_mean_diffCov_', signature, '.png' , sep='' ), total_mm_mean, ylims=ylims_mean, diff=TRUE )

# plot the medians
# total_mm_median <- measure_extractor( list_total_mm, "Median", diff=TRUE )
total_mm_median <- measure_extractor( list_mm, "Median", diff=TRUE )

metagene_measure_builder( paste( list_out_dir[[ 'metagene_median' ]], '/metagene_median_diffCov_', signature, '.png' , sep='' ), total_mm_median, ylims=ylims_median, diff=TRUE )
