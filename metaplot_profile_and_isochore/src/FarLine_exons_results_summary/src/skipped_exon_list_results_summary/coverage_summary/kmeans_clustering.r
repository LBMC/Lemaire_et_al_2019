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
source( paste( kmeans_clustering_dir, 'kmeans_heatmap.r', sep='/' ) )
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
message( paste( 'comp_pair', comp_pair ) )
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

ymax <- NULL
if ( '--hm-ymax' %in% args ) {
  index <- which( args == '--hm-ymax' )
  ymax <- as.numeric( args[ index + 1 ] )
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

# write the data frame of experiment design
write.table( seqDepth_df, file=paste( out_dir, 'seqDepth.tsv', sep='/' ), sep='\t', quote=FALSE, row.names=FALSE )



# load annotations
annot <- read.table( annot_file, header=FALSE, stringsAsFactors=FALSE )
annot[ , 6 ] <- as.character( annot[ , 6 ] )
annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]
#######


#######
cov_enough_id <- list()
cov_df_nullSig_list <- list()
k_res_list <- list()
for ( repl in rep_levels ) {
  message( repl )
  cov_enough_id[[ repl ]] <- c( 1:nrow( annot ) )

  ## load the coverages
  if ( length( comp_pair ) == 1 ) {
    cond <- comp_pair
    bw_file <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]

    # evaluate the annotations with too weak signal
    # message( 'min_max_signal' )
    cov_df <- stranded_coverages( bw_file, annot, seqDepth_df[ seqDepth_df$condition==comp_pair & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=FALSE )
    cov_enough_id[[ repl ]] <- intersect( cov_enough_id[[ repl ]], which( apply( cov_df, 1, max ) > min_max_signal ) )

    # load the working coverage
    # message( 'coverage' )
    cov_df <- stranded_coverages( bw_file, annot, off_set=off_set, seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=harmonization )

  } else if ( length( comp_pair ) == 2 ) {
    bw_file_list <- list()
    normFactor_list <- list()

    for ( cond in comp_pair ) {
      # message( cond )
      bw_file_list[[ cond ]] <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]
      normFactor_list[[ cond ]] <- seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ]
    }

    # evaluate the annotations with too weak signal
    cov_df <- stranded_diffCov( bw_file_list, annot, off_set=off_set, normFactor_list, comp_pair, harmonization=TRUE )
    cov_enough_id[[ repl ]] <- intersect( cov_enough_id[[ repl ]], which( apply( abs( cov_df ), 1, max ) > -1 ) )

    # load the working coverage
    cov_df <- stranded_diffCov( bw_file_list, annot, off_set=off_set, normFactor_list, comp_pair, harmonization=harmonization )
  }

  # record the good annotation and filter out the annotations with too weak signal
  # message( 'filter' )
  rownames( cov_df ) <- annot[ , 4 ]
  cov_df_nullSig_list[[ repl ]] <- cov_df[ ! seq( 1, nrow( cov_df ) ) %in% cov_enough_id[[ repl ]], ]
  cov_df <- cov_df[ cov_enough_id[[ repl ]], ]


  ## compute the k-means clusters
  # message( 'clustering' )
  k_res_list[[ repl ]] <- kmeansp2( as.matrix( cov_df ), nbr_cluster, 1000, 100 )


  ## build heatmap of the cluster with all the exons
  # filename of the image
  # message( 'image' )
  if ( length( comp_pair ) == 1 ) {
    #~ prefix <- paste( c( comp_pair, ''), collapse='-' )
    sample <- paste( comp_pair, 'n', repl, '_', sep='' )
  } else if ( length( comp_pair ) == 2 ) {
    #~ prefix <- 'diffCov'
    sample <- paste( comp_pair[ 1 ], 'n', repl, '-', comp_pair[ 2 ], 'n', repl, '_', sep='' )
  }
	half_name <- paste( suffix, '.png' , sep='')
  heatmap_name <- paste( out_dir, '/heatmap_kmeansp2_', sample, 'allExons_', half_name, sep='' )

  kmeans_heatmap( heatmap_name, nbr_cluster, k_res, cov_df, cov_df_nullSig, ymax=ymax )

}
#######
# message( 'RDS' )

#######
## Save the data in RDS files
if ( length( comp_pair ) == 1 ) {
  prefix <- comp_pair
} else if ( length( comp_pair ) == 2 ) {
  prefix <- paste( comp_pair, collapse='-' )
}
saveRDS( object=list(
k_res_list=k_res_list,
cov_enough_id=cov_enough_id,
min_max_signal=min_max_signal,
annot_file=annot_file,
nbr_cluster=nbr_cluster,
condition=comp_pair
), file=paste( out_RDS_dir, '/cluster_data_', prefix, '_', suffix, '.RDS', sep='' ) )
#######
