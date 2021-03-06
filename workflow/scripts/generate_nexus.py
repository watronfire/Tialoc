import dendropy
import argparse
import numpy as np
import pandas as pd


def load_alignment( loc, verbose=True ):
    if verbose:
        print( "Loading alignment... ", end="" )

    align = dendropy.DnaCharacterMatrix.get( path=loc, schema="fasta" )

    # For reasons unbenounced to me. Somewhere question marks get replaced with with underscores.
    for taxa in align.taxon_namespace:
        taxa.label = taxa.label.replace( "?", "_" )

    if verbose:
        print( "Done" )
        print( "Loaded alignment with {} sequences".format( len( align.taxon_namespace ) ) )

    return align 

def load_tree( loc, align_ns, verbose=True ):
    if verbose:
        print( "Loading tree... ", end="" )
    
    tree = dendropy.Tree.get( path=loc, schema="newick", preserve_underscores=True )

    if verbose:
        print( "Done" )
        print( "Loaded tree with {} tips".format( len( tree.taxon_namespace ) ) )

    assert len( tree.taxon_namespace ) <= len( align_ns ), "Error, tree cannot contain more tips than sequences in alignment. Please use alignment used to generate tree"
    
    # correct taxon_namespace problems
    for i in tree.taxon_namespace:
        i.label = i.label.replace( " ", "_" )

    return tree

def filter_alignment( align, tree, verbose=True ):
    align_u = np.setdiff1d( [i.label for i in align.taxon_namespace], [i.label for i in tree.taxon_namespace] )

    if verbose:
        print( "{} sequences found in alignment but not tree. Removing... ".format( len( align_u ) ), end="" )

    align_ut = [i for i in align.taxon_namespace if i.label in align_u]

    align.discard_sequences( align_ut )
    align.purge_taxon_namespace()    

    if verbose:
        print( "Done" )

    return align

## This works but BEAUTI cannot parse.
def add_metadata( tree, metadata, fields ):
    for leaf in tree.leaf_node_iter():
        entry = metadata.loc[leaf.taxon.label]
        comment = []
        for f in fields:
            leaf.annotations.add_new( f, entry[f] )

def extract_traits( tree, md_loc, fields, output ):
    metadata = pd.read_csv( md_loc, usecols=["strain"] + fields )

    # TODO: this is a brute fix to a larger problem involving special characters in strain names
    metadata["strain"] = metadata["strain"].apply( lambda x: x.replace( "?", "_" ) )

    tree_labels = [i.label for i in tree.taxon_namespace]

    metadata = metadata.loc[metadata["strain"].isin( tree_labels )]
    metadata = metadata.set_index( "strain" )

    metadata.to_csv( output.replace( ".nexus", "_traits.tsv" ), sep="\t" )

    return metadata

def merge_tree_align( args ):
    a = load_alignment( args.alignment )

    t = load_tree( args.tree, a.taxon_namespace )

    if args.fields is not None:
        md = extract_traits( t, args.metadata, args.fields, args.output )
        add_metadata( t, md, args.fields )

    a = filter_alignment( a, t )

    starting_size = len( t.taxon_namespace )

    tl = dendropy.TreeList( [t] )

    ds = dendropy.DataSet( [tl, a] )
    ds.unify_taxon_namespaces()

    if len( ds.taxon_namespaces[0] ) != starting_size:
        print( "Warning: DataSet has a different number of taxons than input tree ({} vs. {})".format( len( ds.taxon_namespaces[0] ), starting_size ) )
    ds.write( path=args.output, schema="nexus" )
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Combines an input alignment and tree into a nexus file for use with BEAST" )
    
    # Input arguments
    parser.add_argument( "-a", "--alignment", required=True, help="input alignment in fasta format" )
    parser.add_argument( "-t", "--tree", required=True, help="input tree in newick format" )
    parser.add_argument( "-m", "--metadata", required=False, help="metadata file for tips in tree" )
    parser.add_argument( "-f", "--fields", nargs="+", help="metadata fields to append to tree" )
    
    # Output arguments
    parser.add_argument( "-o", "--output", required=True, help="output nexus file" )

    arguments = parser.parse_args()

    merge_tree_align( arguments )
