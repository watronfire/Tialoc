import dendropy
import argparse
import numpy as np


def load_alignment( loc, verbose=True ):
    if verbose:
        print( "Loading alignment... ", end="" )

    align = dendropy.DnaCharacterMatrix.get( path=loc, schema="fasta" )
    if verbose:
        print( "Done" )
        print( "Loaded alignment with {} sequences".format( len( align.taxon_namespace ) ) )

    return align 

def load_tree( loc, align_ns, verbose=True ):
    if verbose:
        print( "Loading tree... ", end="" )
    
    tree = dendropy.Tree.get( path=loc, schema="newick" )

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
    

def merge_tree_align( args ):
    a = load_alignment( args.alignment )

    t = load_tree( args.tree, a.taxon_namespace )

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
    
    # Output arguments
    parser.add_argument( "-o", "--output", required=True, help="output nexus file" )

    arguments = parser.parse_args()

    merge_tree_align( arguments )
