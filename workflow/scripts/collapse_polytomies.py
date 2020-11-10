import argparse
from dendropy import Tree

def collapse_polytomies( input_tree, limit, output ):


    count = 0

    bl = [i.edge_length for i in input_tree.internal_nodes()]
    bl = [i for i in bl if i is not None]
    cn = sum( [i < limit for i in bl] )

    print( "{} of {} nodes estimated to be collapsed ({:.2f}%) with a limit of {}".format( cn, len( input_tree.internal_nodes() ), cn / len( input_tree.internal_nodes() ) * 100, args.limit ) )

    input_tree.collapse_unweighted_edges( threshold=limit )

    diff = len( bl ) - len( input_tree.internal_nodes() )
    print( "{} of {} nodes actually were collapsed ({:.2f}%)".format( diff, len( bl ), diff / len( bl ) * 100 ) )

    input_tree.write( path=output, schema="newick" )


def main( arguments ):
    if arguments.url:
        tree = Tree().get( url=arguments.url, schema="newick" )
    elif arguments.path:
        tree = Tree().get( path=arguments.path, schema="newick" )
    else:
        raise ValueError( "No input argument found. Please specify either --url or --path" )

    collapse_polytomies( tree, arguments.limit, arguments.output )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Collapses a tree based on a minimum branch length" )

    # Initialize optional arguments
    parser.add_argument( "-o", "--output", help="ouput location" )
    parser.add_argument( "-l", "--limit", type=float, help="minimum branch length modifier" )

    group = parser.add_mutually_exclusive_group()
    group.add_argument( "-u", "--url", help="URL of input tree" )
    group.add_argument( "-p", "--path", help="location of input tree" )

    args = parser.parse_args()

    main( args )