#!/usr/bin/Rscript

## extend the splicing site regions for window centered on borders
reg_ext_window <- function ( tab_SS, window_size ) {
  tab_SS[ tab_SS$strand == 1, 'start' ] <- tab_SS[ tab_SS$strand == 1, 'start' ] - floor( window_size / 2.0 )
  tab_SS[ tab_SS$strand == 1, 'end' ] <- tab_SS[ tab_SS$strand == 1, 'end' ] + ( ceiling( window_size / 2.0 ) - 1 )

  tab_SS[ tab_SS$strand == -1, 'start' ] <- tab_SS[ tab_SS$strand == -1, 'start' ] - ( ceiling( window_size / 2.0 ) - 1 )
  tab_SS[ tab_SS$strand == -1, 'end' ] <- tab_SS[ tab_SS$strand == -1, 'end' ] + floor( window_size / 2.0 )

  return( tab_SS )
}
