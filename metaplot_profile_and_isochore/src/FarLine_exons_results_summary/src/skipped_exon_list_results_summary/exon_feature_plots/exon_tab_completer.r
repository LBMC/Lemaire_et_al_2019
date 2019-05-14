#!/usr/bin/Rscript

exon_tab_completer_dir <- dirname(sys.frame(1)$ofile)
source( paste( exon_tab_completer_dir, 'val_retriever.r', sep='/' ) )
source( paste( exon_tab_completer_dir, 'order_conv.r', sep='/' ) )
source( paste( exon_tab_completer_dir, 'flanquing_exon_feat.r', sep='/' ) )
source( paste( exon_tab_completer_dir, 'intron_length.r', sep='/' ) )


exon_tab_completer <- function ( exon_tab, feat_tab ) {
  # exon_tab needs columns: 'join_id', 'exons_flanquants', 'id_gene'
  # feat_tab has to be with all the columns that FasterDB provide and the column 'join_id'

  # complete the skipped exon tables with some features
  exon_feat_col <- c( 'START.END', 'EXON_LENGTH', 'INTRON_LENGTH_BEFORE', 'INTRON_LENGTH_AFTER', 'STRENGTH_ACCEPTOR', 'STRENGTH_DONOR', 'TOTAL_EXON_NUMBER', 'STATE', 'FIRST_EXON_ALTERNATIVE', 'LAST_EXON_ALTERNATIVE', 'EXON_DELETION', 'EXON_SKIPPING', 'ACCEPTOR', 'DONOR', 'RETENTION_BEFORE', 'RETENTION', 'X.GC_BEFORE_ACCEPTOR', 'X.GC_AFTER_ACCEPTOR', 'X.GC_BEFORE_DONOR', 'X.GC_AFTER_DONOR', 'FREE_ENERGY_BEFORE_ACCEPTOR', 'FREE_ENERGY_AFTER_ACCEPTOR', 'FREE_ENERGY_ACCEPTOR', 'FREE_ENERGY_BEFORE_DONOR', 'FREE_ENERGY_AFTER_DONOR', 'FREE_ENERGY_DONOR' )
  exon_tab_complete <- merge( exon_tab, feat_tab, by='join_id' )[ , c( colnames( exon_tab ), exon_feat_col ) ]

  # compute the flanquing intron lengths of the skipping events
  # exon_tab_complete[ , 'INTRON_LENGTH_BEFORE' ] <- intron_length( exon_tab_complete, feat_tab, pos='up' )
  # exon_tab_complete[ , 'INTRON_LENGTH_AFTER' ] <- intron_length( exon_tab_complete, feat_tab, pos='down' )

  # recover features for flanquing exons of the skipping events
  exon_feat_tab <- list(
                        up=flanquing_exon_feat( exon_tab_complete, feat_tab, pos='up' ),
                        down=flanquing_exon_feat( exon_tab_complete, feat_tab, pos='down' )
                       )
    # flanquing exons length
    exon_len_col <- 'EXON_LENGTH'
    exon_tab_complete[ , 'EXON_LENGTH_BEFORE' ] <- exon_feat_tab[[ 'up' ]][ , exon_len_col ]
    exon_tab_complete[ , 'EXON_LENGTH_AFTER' ] <- exon_feat_tab[[ 'down' ]][ , exon_len_col ]

    # force acceptor and donor
  exon_flanq_col <- list(
                        up=c( 'FORCE_ACCEPTOR_BEFORE', 'FORCE_DONOR_BEFORE' ),
                        down=c( 'FORCE_ACCEPTOR_AFTER', 'FORCE_DONOR_AFTER' )
                       )

  exon_tab_complete[ , exon_flanq_col[[ 'up' ]] ] <- exon_feat_tab[[ 'up' ]][ , exon_flanq_col[[ 'up' ]] ]
  exon_tab_complete[ , exon_flanq_col[[ 'down' ]] ] <- exon_feat_tab[[ 'down' ]][ , exon_flanq_col[[ 'down' ]] ]

  return( exon_tab_complete )
}
