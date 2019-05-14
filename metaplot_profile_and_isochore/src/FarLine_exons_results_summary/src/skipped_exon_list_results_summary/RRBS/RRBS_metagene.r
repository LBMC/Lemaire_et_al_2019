#!/usr/bin/Rscript

require( 'ggplot2' )
require( 'reshape2' )

## process arguments
args <- commandArgs(trailingOnly = TRUE)
fract_rds_vec <- unlist( strsplit( args[1], split=',' ) )
out_prefix <- args[2]; dir.create( file.path( dirname( out_prefix ) ), showWarnings=FALSE, recursive=TRUE )

prefix_vec <- do.call( rbind, strsplit( fract_rds_vec, split='/' ) )[ , 2 ]
if ( '--prefixes' %in% args ) {
  index <- which( args == '--prefixes' )
  prefix_vec <- unlist( strsplit( args[ index + 1 ], split=',' ) )
}

comp_pair <- NULL
if ( '--comp-pair' %in% args ) {
  index <- which( args == '--comp-pair' )
  comp_pair <- unlist( strsplit( args[ index + 1 ], split=',' ) )
}

mode_graph <- 'raw' #~ available values: 'raw' and 'diff'
if ( '--mode' %in% args ) {
  index <- which( args == '--mode' )
  mode_graph <- args[ index + 1 ]
}

color_pallette <- NULL
if ( '--color-pallette' %in% args ) {
  index <- which( args == '--color-pallette' )
  color_pallette <- paste0( '#', unlist( strsplit( args[ index + 1 ], split=',' ) ) )
}


##### MAIN #####
## load the tables of 5mC fractions
fraction_obj_list <- lapply( fract_rds_vec,
  function ( xxx ) {
    readRDS( xxx )
  }
)

## extract some values
# extract the fraction_list
fraction_list <- lapply(
  fraction_obj_list,
  function ( fraction_obj ) {
    fraction_obj[[ 'fractions' ]]
  }
)

# extract experimental design and mC rate categories
  # ( have to be the same for all the items of fraction_obj_list )
df_expdesign <- fraction_obj_list[[ 1 ]][[ 'df_expdesign' ]]
mC_rate_vec <- rownames( fraction_obj_list[[ 1 ]][[ 'fractions' ]][[ 1 ]] )

if ( is.null( comp_pair ) ) {
  comp_pair <- levels( df_expdesign$condition )
}


## compute the median over the replicates for each condition
for ( mC_rate in mC_rate_vec ) {
  # mC_rate <- ']87:]'

  fractions_med <- list()
  for ( pref_ind in 1:length( prefix_vec ) ) {
    prefix <- prefix_vec[ pref_ind ]
    fractions_med[[ prefix ]] <- list()

    # treat if mode_graph is raw of differential fractions of methCyt
    if ( mode_graph == 'raw' ) {
      sub_comp_pair <- comp_pair
      names( sub_comp_pair ) <- comp_pair
    } else if ( mode_graph == 'diff' ) {
      sub_comp_pair <- paste( comp_pair[ -1 ], comp_pair[ 1 ], sep='-' )
      names( sub_comp_pair ) <- comp_pair[ -1 ]
    }

    # compute the medians per condition
    for ( condition in names( sub_comp_pair ) ) {
      cond_ind <- which( df_expdesign$condition == condition )
      fractions_mC_rate <- list()
      for ( index in cond_ind ) {
        # compute the differential fraction with the reference condition, for each replicate
        if ( mode_graph == 'diff' ) {
          cond_ref_ind <- which( df_expdesign$condition == comp_pair[ 1 ] & df_expdesign$replicate == df_expdesign$replicate[ index ] )
          fraction_obj_list[[ pref_ind ]][[ 'fractions' ]][[ index ]][ mC_rate, ] <- fraction_obj_list[[ pref_ind ]][[ 'fractions' ]][[ index ]][ mC_rate, ] - fraction_obj_list[[ pref_ind ]][[ 'fractions' ]][[ cond_ref_ind ]][ mC_rate, ]
        }

        # recover the (differential) fractions
        repl <- df_expdesign$replicate[ index ]
        fractions_mC_rate[[ repl ]] <- fraction_obj_list[[ pref_ind ]][[ 'fractions' ]][[ index ]][ mC_rate, ]
      }
      # compute the median over the replicates
      fract_med <- as.data.frame( rbind( apply( do.call( rbind, fractions_mC_rate ), 2, median ) ) )
      rownames( fract_med ) <- NULL
      fractions_med[[ prefix ]][[ sub_comp_pair[ condition ] ]] <- fract_med #~t( data.frame( value=fract_med, position=as.numeric( names( fract_med ) ) ) )
    }
  }


  ## build the metagene
  gg_df <- melt( fractions_med )
  gg_df$variable <- as.numeric( as.character( gg_df$variable ) )
  gg_df$L1 <- factor( gg_df$L1, levels=prefix_vec )
  gg_df$L2 <- factor( gg_df$L2, levels=sub_comp_pair )
  leg_nrow <- ( length( levels( factor( gg_df$L1 ) ) ) + 1 ) %/% 2
  # print( str( gg_df ) )

  fig_aes <- aes( x=variable, y=value, color=L1, linetype=L2 )
  if ( length( levels( gg_df$L2 ) ) == 1 ) {
    # figure <- ggplot( gg_df, aes( x=variable, y=value, color=L1 ) )
    fig_aes$linetype <- NULL
  }
  figure <- ggplot( gg_df, fig_aes )
  figure <- figure + xlab( "position to splice site" ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom") + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))

  # set label and limits for y-axis
  if ( mode_graph == 'raw' ) {
    figure <- figure + ylim( 0, 1 ) + ylab( "fraction of exons" )
  } else if ( mode_graph == 'diff' ) {
    figure <- figure + ylim( -1, 1 ) + ylab( "differential fraction of exons" )
  }

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_colour_manual( values=color_pallette )
  }

  # set the suffix of the name of the figure
  if ( mode_graph == 'raw' ) {
    suffix <- 'fract.png'
  } else if ( mode_graph == 'diff' ) {
    suffix <- 'diffFract.png'
  }
  out_file <- paste( out_prefix, mC_rate, suffix, sep='_' )
  # print( out_file )
  png( out_file )
  print(
    figure + geom_line(size=1.5)
    )
  dev.off()

}
