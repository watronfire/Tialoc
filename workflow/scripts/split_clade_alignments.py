import argparse
import argparse
from subprocess import call
import pandas as pd
import tempfile
import os

def load_clades( lineages, clades ):
    df = pd.read_csv( lineages )
    df["group"] = "?"
    for i in clades:
        df.loc[df["lineage"].str.startswith( i ),"group"] = i
    return df


def split_clades( args ):
    lineages = load_clades( args.lineages, args.clades )

    print( args.clades )

    for i in args.clades:
        with tempfile.NamedTemporaryFile( suffix=".txt", mode="w+", delete=False ) as temp_clade:
            entries = lineages.loc[lineages["group"]==i,"taxon"].to_list()
            entries.append( "CHN/Hubei-Wuhan/MN908947.3/2019-12-26" )
            entries.append( "CHN/Hubei-Wuhan/LR757998.1/2019-12-26" )
            entries = list( set( entries ) )
            temp_clade.write( "\n".join( entries ) )


        command = ["module load seqtk",
                   "seqtk subseq {} {} > {}".format( args.alignment, temp_clade.name, os.path.join( args.output, "clade_{}.fasta".format( i ) ) )]

        call( " && ".join( command ), shell=True )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Splits alignment into clades based on pangolin lineage assignments" )

    # Initialize optional arguments
    parser.add_argument( "-a", "--alignment", help="location" )
    parser.add_argument( "-l", "--lineages", help="location" )
    parser.add_argument( "-c", "--clades", nargs="+", help="modifier" )
    parser.add_argument( "-o", "--output", help="location" )
    parser.add_argument( "-r", "--root", nargs="+", help="sequences to include in everytree as a root" )

    arguments = parser.parse_args()
    split_clades( arguments )
