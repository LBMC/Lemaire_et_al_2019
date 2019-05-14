#!/usr/bin/Rscript

## Function to retrieve the value given after a given option
arg_value <- function ( arg_name, args=commandArgs( trailingOnly=TRUE ), type_fc=as.character, default_value="" ) {
  # arg_name is the option (ex: --psi, -a, ... )
  # args is the vector of arguments of the command line
  # type_fc is the function to apply on the value given to the option (ex: as.character, as.numeric, ...)

  if ( arg_name %in% args ) {
    return( type_fc( args[ which( args == arg_name ) + 1 ] ) )
  } else {
    return( default_value )
  }
}

# message( " Function 'args_value' is added to the environment. " )
