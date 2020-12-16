import argparse
from functools import reduce
from dendropy import Tree
import pandas as pd
import os
from fnmatch import fnmatch
import re

"""
Pipeline steps:
 - Load in metadata
 - Load in catchment tree
 - Load in collapsed nodes
 - Load in closest_in_db.csv
 - Catchment tree: Retain one unique sequence per date per country/state.
 - Catchment tree: Retain only earliest sequence from clades that contain >= 90% sequences from one US state or one country. Perform with a preorder traverse. So nodes can't be selected for inclusion if a parent node has been collapsed.
 - Collapsed nodes: Retrain earliest sequence if node is >= 90% from one US state or country.
 - Return list of braches in tree + closest_in_db + interest sequences. 
"""

pd.set_option('display.max_columns', None)

date_re = re.compile( "[0-9]{4}-[0-9]{2}-[0-9]{2}$" )

def load_metadata( file ):
    use_cols = ["strain", "date", "country", "division", "location", "interest"]
    return_df = pd.read_csv( file, usecols=use_cols )
    return_df = return_df.set_index( "strain" )
    return return_df


def load_collapsed_nodes( file ):
    cn = { "name":[], "node":[], "date":[] }
    with open( file, "r" ) as collapsed_nodes:
        next( collapsed_nodes )
        for line in collapsed_nodes:
            line_split = line.strip().split( "," )
            for node in line_split[2].strip( '[]' ).split( " " ):
                # Sometimes a subtree can be listed and I don't know what to do with that at the moment.
                if "subtree" in node:
                    continue

                cn["name"].append( node )
                cn["node"].append( line_split[0] )
                try:
                    cn["date"].append( date_re.search( node )[0] )
                except TypeError:
                    print( "Unable to get date for leaf {} in node {}".format( node, line_split[0] ) )
                    cn["date"].append( "?" )

    cn = pd.DataFrame( cn )
    cn["country"] = cn["name"].apply( lambda x: x.split( "/" )[0] )
    cn["location"] = cn["name"].apply( lambda x: x.split( "/" )[1] )
    cn[["division","location"]] = cn["location"].str.split( "-", n=1, expand=True )

    print( cn.loc[cn["node"]=="collapsed_1267", "country"].value_counts() )


    # Todo: collapsed nodes do not need to be processed in the presense of the tree, so here might be a good place to
    # filter them. Maybe melt and
    return cn


def add_location_fields( tree, md, fields=None ):
    if fields is None:
        fields = ["focus"]

    for i in tree.taxon_namespace:
        i.keep = False
        try:
            entry = md.loc[i.label]
            for j in fields:
                # This may be unnecessary, but might be nice if other attributes need to be added.
                setattr( i, j, entry[j] )
        except KeyError:
            for j in fields:
                setattr( i, j, False )


def prune_redundant_leaves( tree, limit=0.0001, verbose=True, collapse_interest=True ):
    prune = list()
    redundant_c = 0
    missing_c = 0

    for node in tree.internal_nodes():
        visited_locs = dict()
        for child in node.child_node_iter( lambda c: c.is_leaf() and c.edge_length < limit ):
            if child.taxon.interest == "interest":
                if collapse_interest:
                    loc = child.taxon.location if child.taxon.country == "USA" else child.taxon.division
                else:
                    continue
            # there is an odd edge case where the child is a collapsed node, in that case we just prune it, and deal
            # with the problems later.
            elif not child.taxon.interest:
               prune.append( child )
            else:
                loc = child.taxon.division if child.taxon.country == "USA" else child.taxon.country

            try:
                date = date_re.search( child.taxon.label )[0]
            except TypeError:
                missing_c += 1
                prune.append( child.taxon )
                continue

            if loc in visited_locs:
                if date in visited_locs[loc]:
                    redundant_c += 1
                    prune.append( child.taxon )
                else:
                    visited_locs[loc].append( date )
            else:
                visited_locs[loc] = [date]

    if verbose:
        print( "Pruning {} leaves for being redundant.".format( redundant_c ) )
        print( "Pruning {} leaves for missing information.".format( missing_c ) )

    return tree.extract_tree_without_taxa( prune )


def pick_from( taxons, tree, metric, return_nonpicked=False ):
    if metric in ["earliest", "latest"]:
        sorted_taxons = sorted( taxons, key=lambda x: date_re.search( x.label )[0], reverse=(metric=="latest") )
    elif metric in ["basal", "distal"]:
        sorted_taxons = sorted( taxons, key=lambda x: tree.find_node_for_taxon( x ).distance_from_root(), reverse=(metric=="distal") )

    if return_nonpicked:
        returned = sorted_taxons.pop(0)
        return returned, sorted_taxons
    else:
        return sorted_taxons[0]


def subsample_monophylics( tree, limit=0.9, verbose=True, collapse_interest=False, metric="distal" ):
    to_prune = list()

    for node in tree.preorder_internal_node_iter():
        locs = dict()
        contains_interest = False
        for child in node.leaf_iter():
            if child.taxon.interest == "interest":
                loc = child.taxon.location if child.taxon.country == "USA" else child.taxon.division
                contains_interest = True
            else:
                loc = child.taxon.division if child.taxon.country == "USA" else child.taxon.country

            if loc in locs:
                locs[loc].append( child.taxon )
            else:
                locs[loc] = [child.taxon]

        if not collapse_interest and contains_interest:
            continue

        locations_total = sum( len( i ) for i in locs.values() )
        counts = list( map( lambda x: (x, len( locs[x] ) ), locs ) )
        counts = sorted( counts, key=lambda x: x[1] )
        most_common = counts[-1]

        if hasattr( node.parent_node, "collapse" ):
            parent_status = node.parent_node.collapse
        else:
            parent_status = False

        # Need to get around case where a node only has collapsed nodes as children.
        if not most_common[0]:
            node.mcollapse = False
        elif ( most_common[1] / locations_total >= 0.9 ) and not parent_status:
            node.mcollapse = True

            keep = pick_from( locs[most_common[0]], tree, metric=metric )

            for key in locs.keys():
                for value in locs[key]:
                    if value != keep:
                        to_prune.append( value )
        else:
            node.mcollapse = False

    if verbose:
        print( "Pruning {} leaves for being in monophyletic clades.".format( len( to_prune ) ) )
    return tree.extract_tree_without_taxa( to_prune )


def subsample_polytomies( tree, metric="distal", verbose=True ):
    to_prune = list()
    for node in tree.preorder_internal_node_iter():
        locs = dict()
        count = 0
        for child in node.child_node_iter( lambda c: c.is_leaf() ):
            if child.taxon.interest == "interest":
                loc = "interest"
            else:
                loc = child.taxon.division if child.taxon.country == "USA" else child.taxon.country
            if loc in locs:
                locs[loc].append( child.taxon )
            else:
                locs[loc] = [child.taxon]
            count += 1

        if count > 2:
            node.pcollapse = True
            for children in locs.values():
                keep, prune = pick_from( children, tree, metric=metric, return_nonpicked=True )
                to_prune.extend( prune )
                #print( "\t" + "\n\t".join( map( str, locs[i] ) ) )

    if verbose:
        print( "Pruning {} leaves for being at polytomies.".format( len( to_prune ) ) )
    return tree.extract_tree_without_taxa( to_prune )


def parse_llama_output( args ):
    # Load metadata
    md = map( load_metadata, args.metadata )
    md = reduce( lambda a, b: a.combine_first( b ), md ) #  correct shape (81613, 6)

    # Load in catchment trees
    catchment_trees = [i for i in os.listdir( os.path.join( args.llama, "catchment_trees") ) if fnmatch( i, "localsubtree_*.newick")]
    catchment_trees = map( lambda x: Tree.get( path=os.path.join( args.llama, "catchment_trees", x ), schema="newick", label=x ), catchment_trees )

    # Load in closes_in_db
    closest_in_db = pd.read_csv( os.path.join( args.llama, "closest_in_db.csv" ), usecols=["strain"] )
    closest_in_db = closest_in_db["strain"].to_list()
    keep = md.loc[md["interest"]=="interest"].index.to_list()

    report = list()
    for tree in catchment_trees:
        starting_size = len( tree.taxon_namespace )
        add_location_fields( tree, md, fields=["date", "country", "division", "location", "interest"] )
        tree = prune_redundant_leaves( tree, limit=0.00003, collapse_interest=False, verbose=False )
        redundant_leaves = starting_size - len( tree.leaf_nodes() )
        tree = subsample_monophylics( tree, collapse_interest=False, metric="earliest", verbose=False )
        monophyletic_leaves = starting_size - redundant_leaves - len( tree.leaf_nodes() )
        tree = subsample_polytomies( tree, metric="distal", verbose=False )
        polytomy_leaves = starting_size - redundant_leaves - monophyletic_leaves - len( tree.leaf_nodes() )
        tree.purge_taxon_namespace()

        end_size = len( tree.taxon_namespace )

        report.append( "Tree {}: {} to {} nodes ({:.1f}% reduction).".format( tree.label, starting_size, end_size, ( 1 - end_size / starting_size ) * 100 ) )
        report.append( "\t {} leaves pruned for being redundant.".format( redundant_leaves ) )
        report.append( "\t {} leaves pruned for being in monophyletic clades.".format( monophyletic_leaves ) )
        report.append( "\t {} leaves pruned for being at polytomies.".format( polytomy_leaves ) )
        report.append( "" )

        keep.extend( [i.label for i in tree.taxon_namespace] )

    print( "{} leaves remain in the final trees".format( len( keep ) ) )
    print( "\n".join( report ) )
    print( "{} sequences will be included from closest_in_db.csv".format( len( closest_in_db ) ) )

    keep.extend( closest_in_db )
    keep.append( "CHN/Hubei-Wuhan/MN908947.3/2019-12-26" )
    keep.append( "CHN/Hubei-Wuhan/LR757998.1/2019-12-26" )
    unique_keep = list( set( keep ) )

    with open( args.output, "w" ) as outputfile:
        outputfile.write( "\n".join( unique_keep ) )

    print( "{} unique sequences were written to output".format( len( unique_keep ) ) )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="subsamples sequences from llama output catchment files" )

    # Initialize optional arguments
    parser.add_argument( "-l", "--llama", help="location of llama output" )
    parser.add_argument( "-o", "--output", help="location to save list of subsampled sequences" )
    parser.add_argument( "-m", "--metadata", nargs='+', help="metadata for all sequences present. Can provide list." )
    parser.add_argument( "-p", "--preset", type=int, default=1, help="Specify how conservative the pruning should be." )

    arguments = parser.parse_args()

    parse_llama_output( arguments )
