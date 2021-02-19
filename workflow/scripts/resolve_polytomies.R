library( ape )

resolve_polytomies <- function( input, output ) {
  t <- read.tree( file=input )
  t <- multi2di( t )
  write.tree( t, output )
}

args = commandArgs( trailingOnly=TRUE )

resolve_polytomies( args[1], args[2] )