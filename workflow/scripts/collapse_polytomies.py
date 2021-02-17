import argparse
from dendropy import Tree
import pandas as pd
from subprocess import call
from tempfile import NamedTemporaryFile

def collapse_polytomies( input_tree, limit ):

    count = 0

    bl = [i.edge_length for i in input_tree.internal_nodes()]
    bl = [i for i in bl if i is not None]
    cn = sum( [i < limit for i in bl] )

    print( "{} of {} nodes estimated to be collapsed ({:.2f}%) with a limit of {}".format( cn, len( input_tree.internal_nodes() ), cn / len( input_tree.internal_nodes() ) * 100, args.limit ) )

    input_tree.collapse_unweighted_edges( threshold=limit )

    diff = len( bl ) - len( input_tree.internal_nodes() )
    print( "{} of {} nodes actually were collapsed ({:.2f}%)".format( diff, len( bl ), diff / len( bl ) * 100 ) )

    return input_tree

# PareTree doesn't work at all.
def PareTree_wrapper( tree, pruners, output ):
    # save pruners to tempfile
    prune_file = NamedTemporaryFile( suffix=".txt", mode="w+", delete=False )
    print( prune_file.name )
    prune_file.write( "\n".join( pruners ) )

    # save tree to tempfile
    tree_file = NamedTemporaryFile( suffix=".new", mode="w+", delete=False )
    tree.write( file=tree_file, schema="newick", real_value_format_specifier=".4E" )
    print( tree_file.name )

    # run PareTree
    command = f" java -jar workflow/scripts/PareTree1.0.2.jar -del {prune_file.name} -f {tree_file.name}"

    call( command, shell=True )

    prune_file.close()
    tree_file.close()

def rename_tree( tree, md_loc ):
    def _load_gisaid_metadata( loc ):
        return_df = pd.read_csv( loc, sep="\t" )
        return_df = return_df.set_index( "strain" )
        return return_df["gisaid_epi_isl"].to_dict()
    gisaid_dict = _load_gisaid_metadata( md_loc )
    error = 0
    pruners = []
    for i in tree.leaf_node_iter():
        try:
            i.taxon.label = gisaid_dict[i.taxon.label]
        except KeyError:
            error += 1
            #print( f"Unable to find gidaid_id for {i.taxon.label}" )
            pruners.append( i.taxon.label )
    print( f"Renamed {len(tree) - error} of {len(tree)} tips in tree ({(len(tree) - error)/len(tree) * 100:.2f}%)")

    return tree, pruners

def main( arguments ):
    prune = False

    if arguments.url:
        tree = Tree().get( url=arguments.url, schema="newick", preserve_underscores=True )
    elif arguments.path:
        tree = Tree().get( path=arguments.path, schema="newick", preserve_underscores=True )
    else:
        raise ValueError( "No input argument found. Please specify either --url or --path" )

    if arguments.rename is not None:
        tree, prune_list = rename_tree( tree, arguments.rename )

    tree = collapse_polytomies( tree, arguments.limit )

    if prune:
        PareTree_wrapper( tree, prune_list, arguments.output )
    else:
        tree.write( path=arguments.output, schema="newick" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Collapses a tree based on a minimum branch length" )

    # Initialize optional arguments
    parser.add_argument( "-o", "--output", help="ouput location" )
    parser.add_argument( "-l", "--limit", type=float, help="minimum branch length modifier" )
    parser.add_argument( "-r", "--rename", help="gisaid metadata file if optional renaming is to be performed.")

    group = parser.add_mutually_exclusive_group()
    group.add_argument( "-u", "--url", help="URL of input tree" )
    group.add_argument( "-p", "--path", help="location of input tree" )


    args = parser.parse_args()

    main( args )