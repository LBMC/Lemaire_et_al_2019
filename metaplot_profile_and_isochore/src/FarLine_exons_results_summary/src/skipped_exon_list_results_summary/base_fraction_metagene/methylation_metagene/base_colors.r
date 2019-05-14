#!/usr/bin/Rscript

base.colors <- function ( n, start=c( 0, 0, 0 ), end=c( 1, 1 ,1 ), alpha = 1)
{
  red_vec <- seq( from=start[ 1 ], to=end[ 1 ], length.out=n )
  green_vec <- seq( from=start[ 2 ], to=end[ 2 ], length.out=n )
  blue_vec <- seq( from=start[ 3 ], to=end[ 3 ], length.out=n )

  rgb_vec <- rep( '', n )

  for ( xxx in c( 1:length( red_vec ) ) ) {
    rgb_vec[ xxx ] <- rgb( red_vec[ xxx ], green_vec[ xxx ], blue_vec[ xxx ], alpha=alpha )
  }

  return( rgb_vec )
}
