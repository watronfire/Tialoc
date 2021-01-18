import argparse
from dendropy import Tree
from parse_llama_output import load_metadata, add_location_fields, prune_redundant_leaves, subsample_monophylics


def main( args ):
    md = load_metadata( args.metadata )

    tree = Tree.get( path=args.tree, schema="newick" )
    starting_size = len( tree.taxon_namespace )
    add_location_fields( tree, md, fields=["date", "country", "division", "location", "interest"] )



    tree = prune_redundant_leaves( tree, limit=0.00003, collapse_interest=True, verbose=False )
    redundant_leaves = starting_size - len( tree.leaf_nodes() )
    tree = subsample_monophylics( tree, collapse_interest=True, metric="earliest", verbose=False )
    monophyletic_leaves = starting_size - redundant_leaves - len( tree.leaf_nodes() )
    tree.purge_taxon_namespace()
    end_size = len( tree.taxon_namespace )

    report = list()
    report.append( "Tree reduced from {} to {} nodes ({:.1f}% reduction).".format( starting_size, end_size, (1 - end_size / starting_size) * 100 ) )
    report.append( "\t {} leaves pruned for being redundant.".format( redundant_leaves ) )
    report.append( "\t {} leaves pruned for being in monophyletic clades.".format( monophyletic_leaves ) )
    report.append( "" )

    tree.write( path=args.output, schema="newick" )

if __name__ == "__main__" :
    parser = argparse.ArgumentParser( description="description of program" )

    # Initialize optional arguments
    parser.add_argument( "-t", "--tree", help="tree to be pruned" )
    parser.add_argument( "-m", "--metadata", help="metadata to append to tree" )
    parser.add_argument( "-o", "--output", help="location to save output tree" )

    arguments = parser.parse_args()

    main( arguments )
