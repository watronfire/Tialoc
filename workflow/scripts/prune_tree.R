library( ape )
library( Biostrings )

args = commandArgs( trailingOnly=TRUE )
output = args[4]

# Load data
tree <- read.tree( args[1] )
aln <- readDNAMultipleAlignment( args[2], format="fasta" )
md <- read.csv( args[3], sep="\t" )

# Rename tips in tree and remove tips which don't have an entry in tree
tree$tip.label[which(tree$tip.label %in% md$gisaid_epi_isl)] <- as.character( md$strain[match(tree$tip.label[which(tree$tip.label %in% md$gisaid_epi_isl)], md$gisaid_epi_isl)] )
print( cat( length( tree$tip.label[startsWith( tree$tip.label, "EPI")] ), "tips removed because of lack of metadata entry\n" ) )
tree <- drop.tip( tree, tree$tip.label[startsWith( tree$tip.label, "EPI")])

# Remove leaves that aren't in alignment. Output global tree
print( cat( length( setdiff( tree$tip.label, rownames( aln ) ) ), "tips removed because of lack of entry in alignment\n" ) )
tree <- drop.tip( tree, setdiff( tree$tip.label, rownames( aln ) ) )
write.tree( tree, file=file.path( output, "global.tree" ) )

print( cat( length( tree$tip.label ), "tips remaining in tree\n" ) )

# Write tips in tree to file
write( tree$tip.label, file.path( output, "alignment.txt" ) )
#rowmask( aln, invert=TRUE ) <- which( rownames( aln ) %in% tree$tip.label)
#writeXStringSet( x=as( aln, "XStringSet"), filepath=file.path( output, "alignment.fasta" ), format="fasta" )
#rowmask( aln ) <- NULL

# Filter MD to tips in tree.output and write to file
write.csv( subset( md, md$strain %in% tree$tip.label), file=file.path( output, "metadata.csv" ), row.names=FALSE )

# Filter metadata and alignment to query and write to file.
query_md <- subset( md, interest=="interest" )
write.csv( query_md, file=file.path( output, "query.csv" ), row.names=FALSE )
#rowmask( aln, invert=TRUE ) <- which( rownames( aln ) %in% query_md$strain )
#writeXStringSet( x=as( aln, "XStringSet"), filepath=file.path( output, "query.fasta" ), format="fasta" )
write.table( query_md$strain, file.path( output, "query.txt" ), row.names=FALSE, col.names=FALSE, quote=FALSE )
print( cat( nrow( query_md ), "sequences added to query.\n" ) )
