#!/usr/bin/Rscript

## Script to cluster annotation by kmeans++ arlgorithm, build a heatmap according to the cluster and build a barplot of all the exons, for each replicate
# the exons below the threshold of minimum signal are added at the end of the clustering heatmap

#######
## library loading
# Load the libraries and custom functions
require(bigWig, quietly=TRUE, warn.conflicts=FALSE )
require( purrr, quietly=TRUE, warn.conflicts=FALSE )

# set the environment
all_args <- commandArgs(trailingOnly = F)
kmeans_clustering_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
source( paste( kmeans_clustering_dir, 'metagene_coverage_functions.r', sep='/' ) )
source( paste( kmeans_clustering_dir, 'kmeansp2.r', sep='/' ) )
source( paste( kmeans_clustering_dir, 'kmeans_rds.r', sep='/' ) )
source( paste( kmeans_clustering_dir, 'kmeans_measure_retriever.r', sep='/' ) )
source( paste( kmeans_clustering_dir, 'kmeans_cov_enough.r', sep='/' ) )
#######


#######
## Recover the input arguments
harmonization <- FALSE
out_dir <- '.'
suffix='noS'
min_max_signal <- -1
off_set <- 0

args <- commandArgs(trailingOnly = TRUE)
# print( args )
# stop( 'bouh' )
bw_files <- unlist( strsplit( args[ which( args == '-bw' ) + 1 ], split=',' ) ) #~ list of the bw files (~ sample) to consider
conditions <- unlist( strsplit( args[ which( args == '-cond' ) + 1 ], split=',' ) ) #~ indicate for each bw_files
replicates <- unlist( strsplit( args[ which( args == '-rep' ) + 1 ], split=',' ) ) #~ replicate for each bw_file
# seqDepth_file <- args[ which( args == '-seqDepth' ) + 1 ] #~ "~/analyses_SEBASTIEN/data/bam_files/MBD-Seq_MCF7_auboeuf/MBD-Seq_auboeuf_MCF7_sequencingDepthsBySample.tsv"
annot_file <- args[ which( args == '-annot' ) + 1] #~ "~/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/RNA-Seq_results/refined_siPP_exon_list/splicingSites/exons_down_siDNMT3bXsiPP_1kbUpDown5pSS.bed"
comp_pair <- unlist( strsplit( args[ which( args == '-comp_pair' ) + 1 ], split=',' ) )
nbr_cluster <- as.numeric( args[ which( args == '-nbr_cluster' ) + 1] )


if ( '-sign' %in% args ) {
  index <- which( args == '-sign' )
  suffix <- args[ index + 1 ]
}
if ( '-out_dir' %in% args ) {
    index <- which( args == '-out_dir' )
  out_dir <- args[ index + 1 ]
  out_RDS_dir <- paste( out_dir, '/clustering_data', sep='/' )
  dir.create( out_dir, recursive=TRUE, showWarnings=FALSE )
  dir.create( out_RDS_dir, recursive=TRUE, showWarnings=FALSE )
}
if ( '-scaleHarm' %in% args ) {
  harmonization <- TRUE
}
if ( '-min_max_signal' %in% args ) {
  index <- which( args == '-min_max_signal' )
  min_max_signal <- as.numeric( args[ index + 1 ] )
}
if ( '-off_set' %in% args ) {
  index <- which( args == '-off_set' )
  off_set <- as.numeric( args[ index + 1 ] )
}
#######


#######
## set the environment
# load the sequencing depth information
df_expdesign <- data.frame( file=bw_files, condition=conditions, replicate=replicates, row.names=seq( 1, length( bw_files ) ) )
df_expdesign$file <- as.character( df_expdesign$file )
rep_levels <- levels( as.factor( df_expdesign$replicate ) )

vec_seqDepth <- sapply( df_expdesign[[ 'file' ]], bw_seqDepth_retriever )
seqDepth_df <- cbind( df_expdesign, sequencing_depth=vec_seqDepth )
seqDepth_df[ 'normFactor' ] <- seqDepth_df[ ,'sequencing_depth' ]/min( seqDepth_df[ ,'sequencing_depth' ] )

# message( '> Experiment design')
# print( seqDepth_df )


# load annotations
annot <- read.table( annot_file, header=FALSE, stringsAsFactors=FALSE )
annot[ , 6 ] <- as.character( annot[ , 6 ] )
annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]
#######


#######
cov_enough_id <- c( 1:nrow( annot ) )
cov_df_nullSig_list <- list()
k_res_list <- list()
for ( repl in rep_levels ) {
  message( repl )
  ## load the coverages
  # evaluate the annotations with too weak signal
  cov_df <- kmeans_measure_retriever( df_expdesign, seqDepth_df, annot, repl, comp_pair, cov_harm=FALSE, diff_harm=TRUE )
  cov_enough_id <- kmeans_cov_enough( comp_pair, cov_df, min_max_signal )

  # load the working statistics
  cov_df <- kmeans_measure_retriever( df_expdesign, seqDepth_df, annot, repl, comp_pair, cov_harm=harmonization, diff_harm=harmonization )


  ## record the good annotation and filter out the annotations with too weak signal
  # message( 'filter' )
  cov_df_nullSig_list[[ repl ]] <- cov_df[ ! seq( 1, nrow( cov_df ) ) %in% cov_enough_id, ]
  cov_df <- cov_df[ cov_enough_id, ]


  ## compute the k-means clusters
  # message( 'clustering' )
  k_res_list[[ repl ]] <- kmeansp2( as.matrix( cov_df ), nbr_cluster, 1000, 100 )


  ## build heatmap of the cluster with all the exons
  kmeans_heatmap( comp_pair, repl, suffix, out_dir, nbr_cluster, k_res_list, cov_df, cov_df_nullSig_list )
}
#######
# message( 'RDS' )

#######
## Save the data in RDS files
kmeans_rds( out_RDS_dir, suffix, k_res_list, cov_enough_id, min_max_signal, annot_file, nbr_cluster, comp_pair )
#######
