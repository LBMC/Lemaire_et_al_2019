#!/usr/bin/Rscript

base_fraction_acf_plot_dir <- dirname(sys.frame(1)$ofile)
source( paste( base_fraction_acf_plot_dir, 'base_fraction_metagene_fig_name.r', sep='/' ) )
source( paste( base_fraction_acf_plot_dir, 'base_colors.r', sep='/' ) )

## build the curves of metagenes
base_fraction_acf_plot_fc <- function( list_per_meta, vec_annot_name, base_vec, out_dir_pref='autoCorrelation', color_pallette=NULL, skew_measure=NULL, prefix='', entropy=FALSE, struct_info=FALSE, ent_base_vec=NULL, regexpr_fig_name=NULL, presence.only=FALSE ) {
  # compute the name of the figure (and create the output folder)
  metagene_fig_name <- base_fraction_metagene_fig_name_fc( base_vec, skew_measure=skew_measure, prefix=prefix, entropy=entropy, struct_info=struct_info, ent_base_vec=ent_base_vec, regexpr_fig_name=regexpr_fig_name, presence.only=presence.only )
  metagene_dir <- dirname( metagene_fig_name )
  acf_dir <- paste0( metagene_dir, '/', out_dir_pref )
  dir.create( acf_dir, showWarnings = FALSE, recursive = TRUE )

  fig_name <- paste0( acf_dir, '/', basename( metagene_fig_name ) )
  # message( fig_name )

  # compute confidence interval
  conf.level <- 0.95
  max_len_sig <- max( unlist( lapply( list_per_meta, length ) ) )
  # ciline <- qnorm( ( 1 - conf.level ) / 2 ) / sqrt( max_len_sig )

  # reshape the list of metagenes in a data.frame
  list_per_meta <- lapply( list_per_meta, function( xxx ) { data.frame( acf=tail( xxx, n=-1 ), lag=seq( from=1, to=length( xxx ) - 1, by=1 ) ) } )
  df_meta_acf <- melt( list_per_meta, id.vars=c( 'lag', 'acf' ) )
  df_meta_acf$L1 <- factor( df_meta_acf$L1, levels=vec_annot_name )
  leg_nrow <- ( length( levels( factor( df_meta_acf$L1 ) ) ) + 1 ) %/% 2

  # build the plot
  figure <- ggplot( data=df_meta_acf, aes( x=lag, y=acf, color=L1 ) )

  # set colors if defined
  if ( ! is.null( color_pallette ) ) {
    figure <- figure + scale_color_manual( values=color_pallette )
  }

  # figure <- figure + geom_bar( stat = "identity", position = "dodge" )
  figure <- figure + geom_hline( yintercept = 0 )
  # figure <- figure + geom_hline( yintercept=abs( ciline ), linetype='dashed', size=1.5 ) + geom_hline( yintercept=-abs( ciline ), linetype='dashed', size=1.5 )
  figure <- figure + geom_line( size=1.5 )

  png( fig_name ) #, width=1200, height=1200 )
  print(
    figure + xlab( 'lag' ) + ylab( 'autocorrelation' ) + theme_bw( base_size = 16 ) + theme(legend.position="bottom")  + theme(legend.title=element_blank()) + theme( legend.text=element_text( size=20 ) ) + guides(col = guide_legend(nrow = leg_nrow))
  )
  bouh <- dev.off()

}
