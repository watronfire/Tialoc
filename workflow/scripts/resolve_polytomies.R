library( ape )

resolve_polytomies <- function( input, output ) {
  t <- read.tree( file=input )
  t <- multi2di( t )
  write.tree( t, output )
  
}

resolve_polytomies( snakemake@input[["tree"]], snakemake@output[["bi_tree"]] )