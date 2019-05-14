#!/usr/bin/Rscript

require( ggplot2, quietly=TRUE )
require( reshape2, quietly=TRUE )
require( gridExtra, quietly=TRUE )

all_args <- commandArgs(trailingOnly = F)
script_dir <- dirname(sub("--file=","",all_args[grep("--file",all_args)]))
# script_dir <- '/home/sebastien/analyses_SEBASTIEN/analyse_EMT-DNMT3B-PP_MDA-MCF7/src/WGBS_methExtract'
source( paste( script_dir, 'list_tab_reader.r', sep='/' ) )
source( paste( script_dir, 'seqTab2NbrCpG_fc.r', sep='/' ) )
source( paste( script_dir, 'nb_CpG_boxplot_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_mat_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_meas_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_cat_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_meta_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_meta_median_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_repart_plot.r', sep='/' ) )
source( paste( script_dir, 'meth_heatmap_plot.r', sep='/' ) )
source( paste( script_dir, 'signal2Acf_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_acf_plot.r', sep='/' ) )
source( paste( script_dir, 'meth_meta_plot_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_meas_boxplot_fc.r', sep='/' ) )
source( paste( script_dir, 'meth_cat_boxplot_fc.r', sep='/' ) )


#### MAIN ####
args <- commandArgs(trailingOnly = TRUE)
meth_file_vec <- unlist( strsplit( args[ 1 ], split=',' ) )
seq_file_vec <- unlist( strsplit( args[ 2 ], split=',' ) )
name_vec <- unlist( strsplit( args[ 3 ], split=',' ) )

cond_vec <- unlist( strsplit( args[ 4 ], split=',' ) )
rep_vec <- unlist( strsplit( args[ 5 ], split=',' ) )
comp_pair <- unlist( strsplit( args[ 6 ], split=',' ) )

start_pos <- as.numeric( args[ 7 ] )
window <- as.numeric( args[ 8 ] )
min_reads <- args[ 9 ]
out_prefix <- args[ 10 ]

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}

# if diff_mode, compute the difference in methylation between the conditions and the first condition, else compute raw %mCpG for all conditions
diff_mode <- FALSE
scale_range <- c( 0, 1 )
if ( '--diff-mode' %in% args ) {
  diff_mode <- TRUE
  scale_range <- c( -1, 1 )
}


## set experimental design data frame of the bisulfite sequencing experiment
df_expdesign <- data.frame( cond_vec, rep_vec )

## build matrix of methylation file names
meth_file_mat <- matrix( meth_file_vec, nrow=length( name_vec ), dimnames=list( name_vec, sample=paste( cond_vec, rep_vec, sep='_' ) ) )

## load seq tables
list_seq_tab <- list_tab_reader( seq_file_vec, name_vec )

## count the nbr of CpG in each sequence
list_nb_CpG <- lapply( list_seq_tab, seqTab2NbrCpG_fc, len.norm=TRUE )

## compute %mCpG through the samples
all_list_meth_meta <- list()
all_list_meth_meas <- list()
all_list_meth_cat <- list()
for ( cond in comp_pair ) {
  if ( diff_mode & cond == comp_pair[ 1 ] ) {
    message( '--- diff & ref' )
    list_ref_cond_mat <- list()
  } else if ( diff_mode & cond != comp_pair[ 1 ] ) {
    all_list_meth_meta[[ cond ]] <- list()
    message( '--- diff & treat' )
  } else {
    all_list_meth_meta[[ cond ]] <- list()
    message( '--- raw' )
  }

  for ( rep in levels( df_expdesign$rep ) ) {
    message( paste0( cond, '_', rep ) )

    ## do nothing if the replicate do not exist
    if ( is.null( df_expdesign[ df_expdesign$cond == cond & df_expdesign$rep == rep, ] ) ) {
      next
    }

    ## initialize the list of mCpG matrix for the ref condition in diff-mode
    if ( diff_mode & cond == comp_pair[ 1 ] ) {
      list_ref_cond_mat[[ rep ]] <- list()
    }

    ## load the tables of methylations
    list_meth_tab <- list_tab_reader( meth_file_mat[ , paste0( cond, '_', rep ) ], name_vec )

    ## build the matrix of methylation rates
    list_meth_mat_list <- list()
    list_meth_meas_list <- list()
    list_meth_cat_list <- list()
    for ( annot_name in names( list_meth_tab ) ) {
      # print( str( list_seq_tab[[ annot_name ]] ) )
      list_meth_mat_list[[ annot_name ]] <- meth_mat_fc( list_meth_tab[[ annot_name ]], list_seq_tab[[ annot_name ]], min.reads=min_reads, start_pos=start_pos, window=window )
      list_meth_meas_list[[ annot_name ]] <- meth_meas_fc( list_meth_tab[[ annot_name ]], list_seq_tab[[ annot_name ]], min.reads=min_reads )
      list_meth_cat_list[[ annot_name ]] <- meth_cat_fc( list_meth_tab[[ annot_name ]], list_seq_tab[[ annot_name ]], min.reads=min_reads )

      ## treat case of diff_mode
      if ( diff_mode & cond == comp_pair[ 1 ] ) {
        list_ref_cond_mat[[ rep ]][[ annot_name ]] <- list_meth_mat_list[[ annot_name ]]
      } else if ( diff_mode & cond != comp_pair[ 1 ] ) {
        list_meth_mat_list[[ annot_name ]] <- list_meth_mat_list[[ annot_name ]] - list_ref_cond_mat[[ rep ]][[ annot_name ]]
      }
    }

    ## do not build plot for ref cond if in diff mode
    if ( diff_mode & cond == comp_pair[ 1 ] ) {
      next
    } else if ( diff_mode & cond != comp_pair[ 1 ] ) {
      suffix <- paste0( cond, 'n', rep, '-', comp_pair[ 1 ], 'n', rep )
    } else {
      suffix <- paste0( cond, 'n', rep )
    }

    ## build the repartition of methylation for the annotations
    meth_repart_plot_fc( list_meth_mat_list, out_prefix, vec_annot_name=name_vec, suffix=suffix )

    ## build the heatmap of methylation
    meth_heatmap_plot_fc( list_meth_mat_list, out_prefix, vec_annot_name=name_vec, suffix=suffix, scale_range=scale_range )

    ## build the autocorrelation plot
    list_meta_acf <- lapply( list_meth_mat_list, signal2Acf_fc )
    meth_acf_plot_fc( list_meta_acf, out_prefix, vec_annot_name=name_vec, out_dir_pref='autoCorrelation', color_pallette=color_pallette, suffix=suffix )

    ## compute the metaplot of methylation
    list_meth_meta_list <- lapply( list_meth_mat_list, meth_meta_fc )
    all_list_meth_meta[[ cond ]][[ rep ]] <- list_meth_meta_list
    all_list_meth_meas[[ cond ]][[ rep ]] <- list_meth_meas_list
    all_list_meth_cat[[ cond ]][[ rep ]] <- list_meth_cat_list

    # break
  }
  # break
}

# saveRDS( all_list_meth_meta, './bouh.RDS' )

# all_list_meth_meta <- readRDS( './bouh.RDS' )
# print( str( all_list_meth_meta ) )

## build various distribution plots
if ( ! diff_mode ) {
  signature <- paste( comp_pair, collapse='_' )

  # nbr of CpG in sequence
  nb_CpG_boxplot_fc( list_nb_CpG, out_prefix, vec_annot_name=name_vec, color_pallette=color_pallette )

  # mean %CpG
  meth_meas_boxplot_fc( all_list_meth_meas, paste0( out_prefix, '_', signature ), vec_annot_name=name_vec, color_pallette=color_pallette )

  # %CpG categories
  meth_cat_boxplot_fc( all_list_meth_cat, paste0( out_prefix, '_', signature ), vec_annot_name=name_vec, color_pallette=color_pallette )
}


## build the figures of the metaplot
all_list_meth_meta_median <- meth_meta_median_fc( all_list_meth_meta )
all_list_meth_meta_median <- modify_depth( .x=all_list_meth_meta_median, .depth=2, .f=function( xxx ) { new <- xxx; colnames( new ) <- c( 'rel_pos', 'value', 'CpG_counts' ); return( new ) } )

for ( cond in comp_pair ) {
  if ( diff_mode & cond == comp_pair[ 1 ] ) {
    next
  } else if ( diff_mode ) {
    signature <- paste0( cond, '-', comp_pair[ 1 ] )
  } else {
    signature <- cond
  }

  meth_meta_plot_fc( all_list_meth_meta_median[[ cond ]], paste0( out_prefix, '_', signature ), vec_annot_name=name_vec, color_pallette=color_pallette, diff.mode=diff_mode )
  # break
}
