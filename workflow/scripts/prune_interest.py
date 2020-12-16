import argparse
from dendropy import Tree
from workflow.scripts.parse_llama_output import load_metadata, add_location_fields, prune_redundant_leaves


def main( args ):
    md = load_metadata( args.metadata )

    tree = Tree.get( path=args.tree, schema="newick" )
    starting_size = len( tree.taxon_namespace )
    add_location_fields( tree, md, fields=["date", "country", "division", "location", "interest"] )

    tree = prune_redundant_leaves( tree, limit=0.00003, collapse_interest=True, verbose=True )

    print( "Pruned tree of {} leaves to {} leaves".format( starting_size, len( tree.leaf_nodes() ) ) )

    tree.write( path=args.output, schema="newick" )

if __name__ == "__main__" :
    parser = argparse.ArgumentParser( description="description of program" )

    # Initialize optional arguments
    parser.add_argument( "-t", "--tree", help="tree to be pruned" )
    parser.add_argument( "-m", "--metadata", help="metadata to append to tree" )
    parser.add_argument( "-o", "--output", help="location to save output tree" )

    arguments = parser.parse_args()

    main( arguments )