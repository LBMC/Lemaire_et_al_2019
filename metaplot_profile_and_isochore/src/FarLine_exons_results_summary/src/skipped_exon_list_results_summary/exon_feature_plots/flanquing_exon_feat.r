#!/usr/bin/Rscript

flanquing_exon_feat_dir <- dirname(sys.frame(1)$ofile)
source( paste( flanquing_exon_feat_dir, 'val_retriever.r', sep='/' ) )
source( paste( flanquing_exon_feat_dir, 'order_conv.r', sep='/' ) )

# extract feature data for flanquing exons
flanquing_exon_feat <- function ( exon_tab_complete, feat_tab, pos ) {
  # exon_tab_complete needs the columns: 'id_gene', 'exons_flanquants'
  # feat_tab needs the columns: 'join_id'
  if ( pos == 'up' ) { field <- 1 } else if ( pos == 'down' ) { field <- 2 }

  # recover exon_gene_pos of flanquing exons
  flanquing_exon_ids <- paste( exon_tab_complete$id_gene, val_retriever( exon_tab_complete[ , 'exons_flanquants' ], sep=':', field=field ), sep='_' )

  # look for exon_gene_pos not related to existing exon
  lack_logi <- ! flanquing_exon_ids %in% feat_tab$join_id
  if ( any( lack_logi ) ) {
    lack_tab <- as.data.frame( matrix( NA, nrow=sum( lack_logi ), ncol=ncol( feat_tab ) ) )
    colnames( lack_tab ) <- colnames( feat_tab )
    lack_tab$join_id <- flanquing_exon_ids[ lack_logi ]
    feat_tab <- rbind( feat_tab, lack_tab )
  }

  # return data of the flanquing exons only
  out_res <- merge( data.frame( join_id=flanquing_exon_ids, stringsAsFactors=F ), feat_tab, by='join_id', all.x=TRUE, sort=FALSE )
  return( out_res )

}
