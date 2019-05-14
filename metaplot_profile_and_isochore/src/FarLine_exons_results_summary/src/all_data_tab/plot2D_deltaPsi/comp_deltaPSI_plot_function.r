#!/usr/bin/Rscript

  # function to retrieve the common and specific rownames of two data.frames
  repart_id <- function ( first_filt, second_filt ) {
    commonID <- rownames( first_filt )[ rownames( first_filt ) %in% rownames( second_filt ) ]
    first_ID <- rownames( first_filt )[ ! rownames( first_filt ) %in% commonID ]
    second_ID <- rownames( second_filt )[ ! rownames( second_filt ) %in% commonID ]
    return( list( commonID=commonID, first_ID=first_ID, second_ID=second_ID ) )
  }

  # function to determine if the name of the pvalue column (glm or Fisher?)
  pval_col_fc <- function ( tableau ) {
    pval_col <- 'pvalue_corrigee_glm'
    if ( ! ( pval_col %in% colnames( tableau ) ) ) {
      pval_col <- 'pvalue_corrigee_Fisher'
    }
    return( pval_col )
  }

	# create canva of the plot
comp_deltaPSI_plot <- function ( first_raw, second_raw, psi_thresh=0.1, pvalue_thresh=0.05, first_name='first', second_name='second', tendency='all', common_only=FALSE, filt_first=TRUE, filt_second=TRUE, pval_highl, correlation=TRUE, enrichment=TRUE, comptage=TRUE ) {
  ## function details
  # "first_raw" and "second_raw" are two data.frames with rownames as the ids of each exon, and containing the column "deltaPSI", "gene_symbol" (the name of the gene containing each exon), "exon_skipped" ( the index of the exon in the containing gene ), and the column of p-values as the last one
  # "psi_thresh" is the lower limit of PSI measure, "pvalue_thresh" is the upper limit of statistical p-value.
  # "first_name" and "second_name" are the name of the two conditions compared in this plot.


  # keep only common exons, and sort the table for siPP and siFactor
  count_exons <- table( c( rownames( first_raw ), rownames( second_raw ) ) )
  shared_exons <- names( count_exons )[ count_exons == 2 ]
  first_tab <- first_raw[ shared_exons, ]
  second_tab <- second_raw[ shared_exons, ]


  # extract the significant exons according to thresholds of deltaPSI and p-value
  #~ first_filt <- first_tab[ !is.na( abs( first_tab$deltaPSI ) ) & abs( first_tab$deltaPSI ) >= psi_thresh & first_tab[ , length( first_tab ) ] <= pvalue_thresh, ]
  #~ second_filt <- second_tab[ !is.na( abs( second_tab$deltaPSI ) ) & abs( second_tab$deltaPSI ) >= psi_thresh & second_tab[ , length( second_tab ) ] <= pvalue_thresh, ]
  first_filt <- first_tab[ !is.na( abs( first_tab$deltaPSI ) ) & abs( first_tab$deltaPSI ) >= psi_thresh, ]
  if ( filt_first ) {
    first_filt <- first_filt[ first_filt[ , pval_col_fc( first_filt ) ] <= pvalue_thresh, ]
  }
  second_filt <- second_tab[ !is.na( abs( second_tab$deltaPSI ) ) & abs( second_tab$deltaPSI ) >= psi_thresh, ]
  if ( filt_second ) {
    second_filt <- second_filt[ second_filt[ , pval_col_fc(second_filt) ] <= pvalue_thresh, ]
  }

  # compute the common ( all, correlated and uncorrelated for deltaPSI ), siPP specific and siFactor specific exons, according to the previous selection of significant exons
  id_list <- repart_id( first_filt, second_filt )

  # select exon on correlation of the deltaPSI sign
  if ( tendency == 'corr' ) {
    id_list <- lapply( id_list, function ( xxx, first_tab, second_tab ) {
      tend_ids <- xxx[ first_tab[ xxx, 'deltaPSI' ] / second_tab[ xxx, 'deltaPSI' ] >= 0 ]
      tend_ids <- tend_ids[ !is.na( tend_ids ) ]
      return( tend_ids )
    }, first_tab=first_tab, second_tab=second_tab )
  } else if ( tendency == 'antiC' ) {
    id_list <- lapply( id_list, function ( xxx, first_tab, second_tab ) {
      tend_ids <- xxx[ first_tab[ xxx, 'deltaPSI' ] / second_tab[ xxx, 'deltaPSI' ] <= 0 ]
      tend_ids <- tend_ids[ !is.na( tend_ids ) ]
      return( tend_ids )
    }, first_tab=first_tab, second_tab=second_tab )

  }

  # print( str( id_list ) )

  # remove the specific exons
  if ( common_only ) {
    id_list[["first_ID"]] <- c()
    id_list[["second_ID"]] <- c()
  }

  # select commonID with good p-value in first table
  if ( pval_highl ) {
    id_list[[ "sig_ID" ]] <- id_list[[ "commonID" ]][ first_tab[ id_list[[ "commonID" ]], pval_col_fc( first_tab ) ] <= pvalue_thresh ]
  }

  commonID <- id_list[["commonID"]]
  first_ID <- id_list[["first_ID"]]
  second_ID <- id_list[["second_ID"]]
  sig_ID <- id_list[["sig_ID"]]

  # compute pearson correlation between siPP and siFactor samples on commonID
  if ( correlation ) {
    cor_error <- tryCatch({
      cor_res <- cor.test( first_filt[ commonID, "deltaPSI" ], second_filt[ commonID, "deltaPSI" ], method='pearson' )
      #~ print( cor_res, stderr() )
      0},
      error=function(e) { message( "correlation test impossible" ); "1"}
      )
  }

  # compute the p-value of enrichment (Fisher's test)
  if ( comptage ) {
    count_table <- matrix( c( length( shared_exons ) - sum( length( first_ID ), length( second_ID ), length( commonID ) ),
    length( first_ID ),
    length( second_ID ),
    length( commonID )
    ), nrow=2, ncol=2,
    dimnames=list( c( paste( "first unsig" ), paste( "first sig" ) ), c( paste( "second unsig" ), paste( "second sig" ) ) ) )
  }

  # print( count_table )
  if ( enrichment ) {
    fisher_error <- tryCatch({
      fisher_res_more <- fisher.test( x=count_table, alternative='greater' )
      fisher_res_less <- fisher.test( x=count_table, alternative='less' )
      #~ print( fisher_res, stderr() )
      0},
      error=function(e) { message( "Fisher test impossible" ); "1"}
      )
  }


  pch=c( 3, 3, 16, 16 )
  colours=c( "#009900", "#0080FF", "#000000", "#990000" )
  cex.xxx <- 1.75
  par( cex=1.25, cex.axis=cex.xxx, cex.lab=cex.xxx, cex.main=cex.xxx, mar=c( 5, 5, 5, 5 ) + 0.1 )

  # plot part
  plot(
        c( 0, 1 ),
        c( 0, 1 ),
        type='n',
        main=paste( '2D plot of deltaPSI in ', first_name, ' and\nin ', second_name, 'for significative exons' ),
        xlab=paste( first_name, 'deltaPSI' ),
        ylab=paste( second_name, 'deltaPSI' ),
        xlim=c( -1, 1 ),
        ylim=c( -1, 1 )
      )
  grid()
  abline( h=0 )
  abline( v=0 )
  abline( a=0, b=1, col="red" )
    # add points with color according to class owner
      # first specific
    points(
          first_tab[  first_ID, "deltaPSI" ],
          second_tab[  first_ID, "deltaPSI" ],
          pch=pch[ 1 ],
          col=colours[ 1 ]
          )
      # second specific
    points(
          first_tab[ second_ID, "deltaPSI" ],
          second_tab[ second_ID, "deltaPSI" ],
          pch=pch[ 2 ],
          col=colours[ 2 ]
          )
      # common
    points(
        first_tab[ commonID, "deltaPSI" ],
        second_tab[ commonID, "deltaPSI" ],
        pch=pch[ 3 ],
        col=colours[ 3 ]
        )
      # significant common
    points(
        first_tab[ sig_ID, "deltaPSI" ],
        second_tab[ sig_ID, "deltaPSI" ],
        pch=pch[ 3 ],
        col=colours[ 4 ]
        )
    # textid <- first_ID
    # text( first_tab[ textid, "deltaPSI" ], second_tab[ textid, "deltaPSI" ], textid, cex=1, pos=2 )

        # add a legend
        text_leg <- c()
        col_leg <- tail( colours, n=-2 )
        pch_leg <- tail( pch, n=-2 )
        if ( ! common_only ) {
          text_leg <- c( text_leg, paste( first_name, "specific" ), paste( second_name, "specific" ) )
          col_leg <- colours
          pch_leg <- pch
        }
        text_leg <- c( text_leg, "common" )
        if ( pval_highl ) {
          text_leg <- c( text_leg, "common sig" )
        }
    legend( -1, 1,
            text_leg,
            pch=pch_leg,
            col=col_leg,
            cex=1.25,
            xjust=0,
            bty='n',
            horiz=TRUE,
            x.intersp=0.4
          )
  #~ dev.off()


    # write the number of events
  #~ count_text <- paste( 'common analysed: ', length( shared_exons ),'\nsig ', first_name, ': ', length( first_ID ), '\ncommon sig: ', length( commonID ), sep='' )
  if ( comptage ) {
    count_text <- ''
    if ( ! common_only ) {
      count_text <- paste( length( shared_exons ),' analysed\n', length( first_ID ) + length( commonID ), ' all sig green\n', length( second_ID ) + length( commonID ), ' all sig blue', sep='' )
    }
    if ( pval_highl ) {
      count_text <- paste( count_text, '\n', length( commonID ), ' common\n', length( sig_ID ), ' common sig', sep='' )
    } else {
      count_text <- paste( count_text, '\n', length( commonID ), ' common sig', sep='' )
    }

    if ( tendency == 'all' ) {
      text( 1, 0, count_text, cex=1.25, pos=2 )
    } else if ( tendency == 'corr' ) {
      text( 0.5, -1, count_text, cex=1.25, pos=3 )
    } else if ( tendency == 'antiC' ) {
      text( -0.5, -1, count_text, cex=1.25, pos=3 )
    }
  }

    # write the pearson correlation results
  if ( correlation ) {
    if ( cor_error == 0 ) {
      text( 0, 1, paste( "pearson correlation:\n-cor: ", round( cor_res$estimate, 2 ), "\n-p-val: ", round( cor_res$p.value, 4 ) ), cex=1, pos=1 )
      trend_line <- list( 'slope'=cov( first_tab[ commonID, "deltaPSI" ], second_tab[ commonID, "deltaPSI" ] )/var( first_tab[ commonID, "deltaPSI" ] ) )
      trend_line[[ 'intercept' ]] <- mean( second_tab[ commonID, "deltaPSI" ] ) - trend_line[[ 'slope' ]] * mean( first_tab[ commonID, "deltaPSI" ] )
      #~ abline( b=cor_res$estimate, a=0, lty=2 )
      abline( b=trend_line[[ 'slope' ]], a=trend_line[[ 'intercept' ]], lty=2 )
    }
  }

    # write the significance of enrichment test
  if ( enrichment ) {
    if ( fisher_error == 0 ) {
      text( 0, -1, paste( "enrichment (Fisher)\nmore common pval: ", round( fisher_res_more$p.value, 4 ), '\nmore spec pval: ', round( fisher_res_less$p.value, 4 ) ), cex=1, pos=3 )
    }
  }

}

#~_dev.off()
