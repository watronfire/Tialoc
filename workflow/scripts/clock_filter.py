import argparse
from datetime import datetime
from scipy import stats
import numpy as np
import argparse
import re
import pandas as pd
from dendropy import Tree

def date_to_dec( date_str ):
    date_obj = datetime( *map( int, date_str.split( "-" ) ) )
    return ( float( date_obj.strftime( "%j" ) ) - 1 ) / 365.25 + float( date_obj.strftime( "%Y" ) )

def prune_clock_outliers( tree, iqd=3, clock_rate=None ):
    date_re = re.compile( "[0-9]{4}-[0-9]{2}-[0-9]{2}$" )

    prune = { "name" : [], "distance" : [], "date" : [] }

    to_prune = []

    for leaf in tree.leaf_node_iter():
        try:
            date = date_re.search( leaf.taxon.label )[0]
        except TypeError:
            to_prune.append( leaf.taxon.label )
            continue
        prune["name"].append( leaf.taxon.label )
        prune["distance"].append( leaf.distance_from_root() )
        prune["date"].append( date_re.search( leaf.taxon.label )[0] )
    prune = pd.DataFrame( prune )
    prune["date"] = prune["date"].apply( date_to_dec )

    # remove bullshit
    to_prune.extend( prune.loc[prune["date"]<2019.0, "name"].to_list() )
    prune = prune.loc[prune["date"]>2019]

    # Compute regression line
    if clock_rate is None:
        print( "Clock rate not specified, infering from data. Rerooting will not be performed." )
        slope, intercept, r_value, p_value, std_err = stats.linregress(x=prune["date"],y=prune["distance"])
        print( "Calculated substitution rate: {}".format( slope ) )
    else:
        slope = clock_rate
        intercept = ( prune["distance"] - prune["date"] * clock_rate ).mean()

    prune["residual"] = prune["distance"] - ( prune["date"] * slope + intercept )
    q3, q1 = np.percentile(prune["residual"], [75 ,25])
    interquartile = q3 - q1
    prune["include"] = abs( prune["residual"] ) < interquartile * iqd

    to_prune.extend( prune.loc[~prune["include"],"name"].to_list() )
    print( "Pruning {} leaves violating clock rate.".format( len( to_prune ) ) )

    return tree.extract_tree_without_taxa_labels( to_prune )


def main( args ):

    tree = Tree.get( path=args.input, schema="newick", preserve_underscores=True )

    tree = prune_clock_outliers( tree, iqd=args.iqd, clock_rate=args.clock_rate )

    tree.write( path=args.output, schema="newick" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Prunes large tree to focal sequences and ancestors" )

    # Initialize optional arguments
    parser.add_argument( "-i", "--input", help="tree to be pruned" )
    parser.add_argument( "-o", "--output", help="location to save output tree" )
    parser.add_argument( "-c", "--clock-rate", type=float, default=None, help="Clock rate to use for clock-based filter. If not specified substitution rate will be infered." )
    parser.add_argument( "--iqd", type=int, default=3, help="Sequences with residual greater than this many iqds will be filtered" )

    arguments = parser.parse_args()

    main( arguments )