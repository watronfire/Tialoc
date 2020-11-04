import argparse
import os
import pandas as pd
import numpy as np
from Bio import AlignIO
from dendropy import Tree

def main( args ):
    # Load metadata and generate strain - gisaid id dictionary.
    md = pd.read_csv( args.metadata, sep="\t" )
    print( "Entries in metadata: {}".format( len( md ) ) )
    md_gisaid_list = md["gisaid_epi_isl"].to_list()
    gisaid_dict = md.loc[~md["gisaid_epi_isl"].isna(), ["strain", "gisaid_epi_isl"]]
    gisaid_dict = gisaid_dict.set_index( "gisaid_epi_isl" )
    gisaid_dict = gisaid_dict["strain"].to_dict()

    # Load alignment
    alignment = AlignIO.read( args.alignment, "fasta" )
    print( "Sequences in alignment: {}".format( len( alignment ) ) )
    alignment_list = [i.name for i in alignment]

    # Load tree
    tree = Tree.get( path=args.tree, schema="newick" )
    print( "Leaves in tree: {}".format( len( tree.taxon_namespace ) ) )

    # Determine leaves which names cannot be assigned.
    tree_leaves = [i.label for i in tree.taxon_namespace]
    tree_leaves = [i.replace( " ", "_" ) for i in tree_leaves]
    leaf_missing_md = np.setdiff1d( tree_leaves, md_gisaid_list )

    # Remove leaves identified
    print( "Leaves in tree but not in metadata: {}".format( len( leaf_missing_md ) ) )
    tree = tree.extract_tree_without_taxa_labels( [i.replace( "_", " " ) for i in leaf_missing_md] )
    tree.purge_taxon_namespace()
    print( "Leaves in tree after pruning: {}".format( len( tree.taxon_namespace ) ) )

    # Rename leaves to match metadata and alignment
    print( "Renaming leaves to match metadata and alignment... ", end="" )
    leaves = list()
    for i in tree.taxon_namespace:
        i.label = gisaid_dict[i.label.replace( " ", "_" )]
        leaves.append( i.label )
    print( "Done" )

    # Remove leaves that aren't in alignment
    leaf_missing_align = np.setdiff1d( leaves, alignment_list )
    print( "Leaves in tree but not in alignment: {}".format( len( leaf_missing_align ) ) )
    tree = tree.extract_tree_without_taxa_labels( leaf_missing_align )
    tree.purge_taxon_namespace()
    print( "Leaves in tree after pruning: {}".format( len( tree.taxon_namespace ) ) )

    # Update tree_leaves list
    tree_leaves = [i.label for i in tree.taxon_namespace]

    # Filter alignment to tips in tree
    tree_alignment = list()
    for i in alignment:
        if i.name in tree_leaves:
            tree_alignment.append( i )
    tree_alignment = AlignIO.MultipleSeqAlignment( tree_alignment )

    # Filter metadata to tips in tree
    tree_md = md.loc[md["strain"].isin( [i.name for i in tree_alignment] )]

    # Filter metadata and alignment to query
    query_md = md.loc[md["interest"] == "interest"]
    interests = query_md["strain"].to_list()
    query_alignment = [i for i in alignment if i.name in interests]
    query_alignment = AlignIO.MultipleSeqAlignment( query_alignment )

    # Write files to disk
    tree.write( path=os.path.join( args.outdir, "global.tree" ), schema="newick" )
    AlignIO.write( tree_alignment, os.path.join( args.outdir, "alignment.fasta" ), "fasta" )
    tree_md.to_csv( os.path.join( args.outdir, "metadata.csv" ), index=False )
    AlignIO.write( query_alignment, os.path.join( args.outdir, "query.fasta" ), "fasta" )
    query_md.to_csv( os.path.join( args.outdir, "query.csv" ), index=False )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Prunes input tree to match alignment, metadata, and focus" )

    # Initialize optional arguments
    parser.add_argument( "-o", "--outdir", help="ouput directory" )
    parser.add_argument( "-m", "--metadata", help="input metdata containing at least name, gisaid_id, and interest field" )
    parser.add_argument( "-t", "--tree", help="input tree containing GISAID ids as tip names" )
    parser.add_argument( "-a", "--alignment", help="input filtered and masked alignment" )

    arguments = parser.parse_args()

    main( arguments )
