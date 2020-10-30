import argparse
from dendropy import Tree

def collapse_polytomies( input_tree, limit, output ):
    tree = Tree().get( url=input_tree, schema="newick" )

    count = 0

    bl = [i.edge_length for i in tree.internal_nodes()]
    cn = sum( [i < limit for i in bl] )

    print( "{} of {} nodes estimated to be collapsed ({:.2f}%)".format( cn, len( tree.internal_nodes() ), cn / len( tree.internal_nodes() ) * 100 ) )

    tree.collapse_unweighted_edges( threshold=limit )

    diff = len( bl ) - len( tree.internal_nodes() )
    print( "{} of {} nodes actually were collapsed ({:.2f}%)".format( diff, len( bl ), diff / len( bl ) * 100 ) )

    tree.write( path=output, schema="newick" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Collapses a tree based on a minimum branch length" )

    # Initialize positional arguments
    parser.add_argument( "input_tree", help="URL of tree to be pruned" )

    # Initialize optional arguments
    parser.add_argument( "-o", "--output", help="ouput location" )
    parser.add_argument( "-l", "--limit", type=float, help="minimum branch length modifier" )

    args = parser.parse_args()

    collapse_polytomies( args.input_tree, args.limit, args.output )