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
## Recover the input arguments
harmonization <- FALSE
out_dir <- '.'
signature <- 'noS'
off_set <- 0

args <- commandArgs(trailingOnly = TRUE)
# print( args )
bw_files <- unlist( strsplit( args[ which( args == '-bw' ) + 1 ], split=',' ) ) #~ list of the bw files (~ sample) to consider
conditions <- unlist( strsplit( args[ which( args == '-cond' ) + 1 ], split=',' ) ) #~ indicate for each bw_files
replicates <- unlist( strsplit( args[ which( args == '-rep' ) + 1 ], split=',' ) ) #~ replicate for each bw_file

annot_file_vec <- unlist( strsplit( args[ which( args == '-annot' ) + 1 ], split=',' ) ) #~ bed file of the annotations to consider
prefixes <- unlist( strsplit( args[ which( args == '-prefixes' ) + 1 ], split=',' ) )
comp_pair <- unlist( strsplit( args[ which( args == '-comp_pair' ) + 1 ], split=',' ) )
message( paste( 'comp_pair', comp_pair ) )


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
  # ylims_list[[ 'mean' ]] <- c( NA, as.numeric( args[ index + 1 ] ) )
  ylims_list[[ 'mean' ]] <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}
if ( '--ylims_median' %in% args ) {
  index <- which( args == '--ylims_median' )
  # ylims_list[[ 'median' ]] <- c( NA, as.numeric( args[ index + 1 ] ) )
  ylims_list[[ 'median' ]] <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}
if ( '--scaleHarm' %in% args ) {
  harmonization <- TRUE
  ylims_list <- list( mean=c( 0, 1 ), median=c( 0, 1 ) )
}


hm_ymax <- NULL
if ( '--hm-ymax' %in% args ) {
  index <- which( args == '--hm-ymax' )
  hm_ymax <- as.numeric( args[ index + 1 ] )
}

ref_add_chr <- FALSE
txt_arg <- '--ref-add-chr'
if ( txt_arg %in% args ) {
  ref_add_chr <- TRUE
}

various_length <- FALSE
nbins <- 0
side_cut <- c( 0, 0 )
ext_up <- 0
ext_dw <- 0
txt_arg <- '--various-length'
if ( txt_arg %in% args ) {
  various_length <- TRUE
  index <- which( args == txt_arg )
  nbins <- as.numeric( args[ index + 1 ] )

  txt_arg <- '--side-cut'
  if ( txt_arg %in% args ) {
    index <- which( args == txt_arg )
    side_cut <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
  }

  txt_arg <- '--ext-up'
  if ( txt_arg %in% args ) {
    index <- which( args == txt_arg )
    ext_up <- as.numeric( args[ index + 1 ] )
  }

  txt_arg <- '--ext-dw'
  if ( txt_arg %in% args ) {
    index <- which( args == txt_arg )
    ext_dw <- as.numeric( args[ index + 1 ] )
  }
}

# differential <- FALSE
# if ( '-diff' %in% args ) {
#   differential <- TRUE
# }

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

gene_concat <- FALSE
if ( '--gene-concat' %in% args ) {
  gene_concat <- TRUE
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
  annot <- annot[ ! grepl( '^#.*', annot[ , 1 ] ), ]
  annot[ , 6 ] <- as.character( annot[ , 6 ] )
  annot <- bedCheckCoords( annot, clean_dup=TRUE, stranded=TRUE )[[ "clean" ]]
  return( annot )
} )
#########


#########
## retrive the coverages and compute the resuming values
fc_vec <- c( 'mean', 'median' )
list_total_mm <- list()
list_total_max <- list()
list_total_mean <- list()
for ( prefix in prefixes ) {
  list_total_max[[ prefix ]] <- list()
  list_total_mean[[ prefix ]] <- list()
  message( paste0( '--- ', prefix ) )
  annot_tab <- annot_list[[ prefix ]]
  for ( cond in comp_pair ) {
    message( paste0( '___ ', cond ) )
    list_mm <- list()
    list_total_max[[ prefix ]][[ cond ]] <- list()
    list_total_mean[[ prefix ]][[ cond ]] <- list()

    for ( repl in df_expdesign$replicate[ df_expdesign$condition == cond ] ) {
      # message( repl )
      bw_file <- df_expdesign$file[ df_expdesign$condition == cond & df_expdesign$replicate == repl ]

      # check refs between BigWig file and annotations
      check_ref_annotVsBw <- function( annot, bw_file, ref_add_chr=F ) {
        bw_reader <- load.bigWig( bw_file )

        annot_ref_col <- annot[ , 1 ]
        if ( ref_add_chr ) {
          annot_ref_col <- paste0( 'chr', annot_ref_col )
        }
        annot_ref <- unique( annot_ref_col )
      
        mis_ref <- setdiff( annot_ref, bw_reader$chroms )
        if ( length( mis_ref ) == length( annot_ref ) ) {
          message( 'All ref missing!' )
          message( paste( annot[ 1, ], collapse='\t' ) )
          stop()
        } else if ( length( mis_ref ) != 0 ) {
          message( 'Several ref missing!' )
          message( paste( mis_ref, collapse='\t' ) )
          message( 'Filter the annotations!' )
          annot <- annot[ which( ! sapply( annot_ref_col, function( xxx, yyy ) { xxx %in% yyy }, mis_ref ) ), ]
        }

      return( annot )
      }

      annot_tab <- check_ref_annotVsBw( annot_list[[ prefix ]], bw_file, ref_add_chr=ref_add_chr )
      

      # retrieve coverage from bam file for all annotations and apply normalization
      if ( various_length ) {
        cov_df <- stranded_coverages_varLen( bw_file, annot_tab, nbins=nbins, off_set=off_set, seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=harmonization, ref_add_chr=ref_add_chr, gene_concat=gene_concat )
        if ( ext_up & ! gene_concat ) {
          annot_out <- annot_tab
          annot_out[ , 3 ] <- annot_out[ , 2 ]
          annot_out[ , 2 ] <- annot_out[ , 2 ] - ext_up
	  log_strand <- annot_out[ , 6 ] %in% c( '-1', '-' )
          annot_out[ log_strand, 2 ] <- annot_out[ log_strand, 3 ]
          annot_out[ log_strand, 3 ] <- annot_out[ log_strand, 2 ] + ext_dw
          cov_out_df <- as.matrix( stranded_coverages( bw_file, annot_out, off_set=0, seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=harmonization, ref_add_chr=ref_add_chr ) )
          cov_df <- cbind( cov_out_df, cov_df )
          cov_out_df <- NULL
          gc()
        }
        # print( str( cov_df ) )
        if ( ext_dw & ! gene_concat ) {
          annot_out <- annot_tab
          annot_out[ , 2 ] <- annot_out[ , 3 ]
          annot_out[ , 3 ] <- annot_out[ , 3 ] + ext_dw
          log_strand <- annot_out[ , 6 ] %in% c( '-1', '-' )
          annot_out[ log_strand, 3 ] <- annot_out[ log_strand, 2 ]
          annot_out[ log_strand, 2 ] <- annot_out[ log_strand, 3 ] - ext_dw
	  cov_out_df <- as.matrix( stranded_coverages( bw_file, annot_out, off_set=0, seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=harmonization, ref_add_chr=ref_add_chr ) )
          cov_df <- cbind( cov_df, cov_out_df )
          cov_out_df <- NULL
          gc()
        }
        # print( str( cov_df ) )
        colnames( cov_df ) <- 1:ncol( cov_df )
        # print( str( cov_df ) )

        if ( side_cut[ 1 ] & ! ext_up ) {
          cov_df <- cov_df[ , ( side_cut[ 1 ] + 1 ):ncol( cov_df ) ]
        }
        if ( side_cut[ 2 ] & ! ext_dw ) {
          cov_df <- cov_df[ , 1:( ncol( cov_df ) - side_cut[ 2 ] ) ]
        }

      } else {
        cov_df <- stranded_coverages( bw_file, annot_tab, off_set=off_set, seqDepth_df[ seqDepth_df$condition==cond & seqDepth_df$replicate==repl, 'normFactor' ], harmonization=harmonization, ref_add_chr=ref_add_chr )
      }

      # repart of the annotations for their maximum coverage
      max_vec <- apply( cov_df, 1, max, na.rm=TRUE )
      list_total_max[[ prefix ]][[ cond ]][[ repl ]] <- max_vec
      filename <- paste( list_out_dir[[ 'repartMax' ]], '/repartMax_', prefix, '_', cond, '_', repl, '_', signature, '.png', sep='' )
      repart_builder( filename, max_vec )


      # repart of the annotations for their average coverage
      mean_vec <- apply( cov_df, 1, mean, na.rm=TRUE )
      list_total_mean[[ prefix ]][[ cond ]][[ repl ]] <- mean_vec
      repart_builder( filename, mean_vec )

      # build the heatmap of coverages
      heatmap_builder( paste( list_out_dir[[ 'heatmapSimple' ]], '/heatmapSimple_', prefix, '_', cond, '_', repl, '_', signature, '.png', sep='' ), cov_df, Rowv=NA, ymax=hm_ymax, sum.order=TRUE )


      # resume the coverages per position
      for ( type_fc in fc_vec ) {
        list_mm[[ type_fc ]][[ repl ]] <- apply( cov_df, 2, eval( parse( text=type_fc ) ) )
      }
      cov_df <- NULL
      gc()
    }
    # message( 'sub' )
    # print( str( list_mm[[ 'mean' ]] ) )

    # resume the replicates the median applied on them for each position
    for ( type_fc in fc_vec ) {
      list_total_mm[[ type_fc ]][[ cond ]][[ prefix ]] <- apply( do.call( rbind, list_mm[[ type_fc ]] ), 2, median )
    }
  }
}

## build repartition plot
filename <- paste( list_out_dir[[ 'repartMax' ]], '/repartMax_', signature, '_', cond, '.png', sep='' )
coverage_repart_plot_fc( list_total_max, vec_annot_name=prefixes, color_pallette=color_pallette, filename=filename )
filename <- paste( list_out_dir[[ 'repartMean' ]], '/repartMean_', signature, '_', cond, '.png', sep='' )
coverage_repart_plot_fc( list_total_mean, vec_annot_name=prefixes, color_pallette=color_pallette, filename=filename )


# structure the resuming value per annotation list in one table
list_total_mm <- list_total_mm %>% modify_depth( .depth=2, .f=function ( xxx ) { do.call( rbind, xxx ) } )
# message( 'total' )
# print( str( list_total_mm ) )
# saveRDS( list_total_mm, 'bouh.RDS' )
#########

# list_total_mm <- readRDS( 'bouh.RDS' )
# message( 'rds' )


#########
## build the figures
for ( type_fc in c( 'mean', 'median' ) ) {
  message( paste0( 'summary type: ', type_fc ) )
    # curve
  total_mm <- list_total_mm[[ type_fc ]] %>% melt()
  colnames( total_mm ) <- c( 'annot', 'X1', 'value', 'condition' )
  total_mm[ , 'annot' ] <- factor( total_mm[ , 'annot' ], levels=prefixes )

  metagene_measure_builder( paste( list_out_dir[[ paste( 'metagene_', type_fc, sep='' ) ]], '/metagene_', type_fc, '_', signature, '.png' , sep='' ), total_mm, ylims=ylims_list[[ type_fc ]], colour=total_mm$annot, color_pallette=color_pallette )

    # heatmap
  if ( length( annot_list ) > 1 ) {
    for ( cond in comp_pair ) {
      heatmap_builder( paste( list_out_dir[[ 'heatmapSimple' ]], '/heatmapSimple_', type_fc, '_', signature, '_', cond, '.png', sep='' ), list_total_mm[[ type_fc ]][[ cond ]], Rowv=NULL )
    }
  }
}
# stop( 'bouh' )
#########
