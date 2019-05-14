#!/usr/bin/Rscript

require( 'ggplot2', quietly = TRUE )
require( 'reshape2', quietly = TRUE )

all_args <- commandArgs(trailingOnly = F)
base_fraction_metagene_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))

## Some function
source( paste( base_fraction_metagene_dir, 'list_tab_reader.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'seqTab2Entropy.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'seqTab2FractMat.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'mat2Metagene.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'listMeta2DF.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'base_fraction_metagene_plot.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'base_fraction_heatmap_plot.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'base_fraction_repart_plot.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'signal2Acf_fc.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'acfOnVec_fc.r', sep='/' ) )
source( paste( base_fraction_metagene_dir, 'base_fraction_acf_plot_fc.r', sep='/' ) )

#####



#### MAIN ####
## process arguments
args <- commandArgs(trailingOnly = TRUE)
# print( paste( args, collapse=' ' ) )

vec_annot_file <- unlist( strsplit( args[1], split=',' ) ) #~ vector of path to the file of sequences (col: chr:start-pos, strand, sequence) [ chr:start-pos and strand used as ids ]
vec_annot_name <- unlist( strsplit( args[2], split=',' ) ) #~ vector of the names for display corresponding to the files of sequences ( submitted in the same order )

window_size <- as.numeric( args[3] ) #~ nbr of bases for the size of the sliding window on the sequence
start_pos <- as.numeric( args[4] ) #~ custom number of the starting position (-700, 0, or 20, ...)
prefix <- args[5] #~ prefix for naming the output figure file

sep_str <- ','
use.regexpr <- FALSE
regexpr_fig_name <- NULL
scale_range <- c( 0, 1 )
if ( '--regexpr' %in% args ) {
  use.regexpr <- TRUE
  sep_str <- '_ppp_'
  index <- which( args == '--regexpr' )
  regexpr_fig_name <- args[ index + 1 ]
}

entropy <- FALSE #~ choose te compute entropy of the DNA bases instead of base fraction or skewness
if ( '--entropy' %in% args ) {
  entropy <- TRUE
  if ( '--entropy' %in% args ) {
    index <- which( args == '--entropy' )
    ent_base_vec <- unlist( strsplit( args[ index + 1 ], split=sep_str ) )
  }


  struct_info <- FALSE #~ choose te compute struct_info of the DNA bases instead of pure entropy
  if ( '--struct-info' %in% args ) {
    struct_info <- TRUE
  }
} else {
  base_vec <- unlist( strsplit( args[6], split=sep_str ) ) #~ c( "G", "C" ); bases for which to compute the fraction in the sliding window

  skew_measure <- NULL #~ bases to consider in addition of base_vec to compute the fraction of bases given in base_vec, default is NULL for basic fraction (0>1), else skewness (-1>1)
  if ( '--skewness' %in% args ) {
    index <- which( args == '--skewness' )
    skew_measure <- unlist( strsplit( args[ index + 1 ], split=sep_str ) )
    scale_range <- c( -1, 1 )
  }
}

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

presence.only <- FALSE
if ( '--presence-only' %in% args ) {
  presence.only <- TRUE
}

ylim_arg <- NULL
if ( '--ylim' %in% args ) {
  index <- which( args == '--ylim' )
  ylim_arg <- as.numeric( unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}


## load the tables of annotations
list_tab <- list_tab_reader( vec_annot_file, vec_annot_name )

if ( entropy ) { # entropy mode
  ## compute the counts
  list_metagene <- lapply( list_tab, seqTab2Entropy, ent_base_vec=ent_base_vec, start_pos=start_pos, window_size=window_size, struct_info=struct_info )

} else { # metagene mode
  ## compute the base fractions
  list_fraction_mat <- lapply( list_tab, seqTab2FractMat, start_pos=start_pos, base_vec=base_vec, window_size=window_size, all.pos=TRUE, skew_measure=skew_measure, use.regexpr=use.regexpr, presence.only=presence.only )
  list_fraction_mat_first_pos <- lapply( list_tab, seqTab2FractMat, start_pos=start_pos, base_vec=base_vec, window_size=window_size, all.pos=FALSE, skew_measure=skew_measure, use.regexpr=use.regexpr, presence.only=presence.only )

  ## compute the metagenes of base fractions
  list_metagene <- lapply( list_fraction_mat, mat2Metagene, reduc_fc=function( xxx ) { mean( xxx, na.rm=TRUE ) } )

  ## compute the metagene of autocorrelation made on each annotation
  list_meta_acf <- lapply( list_fraction_mat_first_pos, signal2Acf_fc )

  ## compute the autocorrelation on the metagenes
  list_acf_on_meta <- lapply( list_metagene, acfOnVec_fc )

  ## build the repartition plots
  base_fraction_repart_plot_fc( list_fraction_mat, vec_annot_name, base_vec, color_pallette=color_pallette, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )

  ## build the heatmaps
  base_fraction_heatmap_plot_fc( list_fraction_mat, vec_annot_name, base_vec, color_pallette=color_pallette, scale_range=scale_range, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )

  ## build meta autocorrelation figures
  base_fraction_acf_plot_fc( list_meta_acf, vec_annot_name, base_vec, out_dir_pref='autoCorrelation', color_pallette=color_pallette, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )

  ## build figures of autocorrelation made on metagenes
  base_fraction_acf_plot_fc( list_acf_on_meta, vec_annot_name, base_vec, out_dir_pref='autoCorOnMeta', color_pallette=color_pallette, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )
}

#### build the figures
## build the curves of metagenes
base_fraction_metagene_plot_fc(
  list_metagene,
  vec_annot_name,
  base_vec,
  color_pallette=color_pallette,
  skew_measure=skew_measure,
  prefix=prefix,
  entropy=entropy,
  struct_info=struct_info,
  ent_base_vec=ent_base_vec,
  regexpr_fig_name=regexpr_fig_name,
  presence.only=presence.only,
  ylim_arg=ylim_arg
)
