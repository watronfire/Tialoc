import argparse
import argparse
from subprocess import call
import pandas as pd
import tempfile
import os

def load_clades( lineages ):
    df = pd.read_csv( lineages )
    
    df.loc[df["lineage"] == "A","group"] = "A"
    df.loc[df["lineage"].str.startswith( "A." ),"group"] = "A"

    df.loc[df["lineage"] == "B","group"] = "B"
    df.loc[df["lineage"].str.startswith( "B." ),"group"] = "B"

    df.loc[df["lineage"]=="B.1","group"] = "B.1"
    df.loc[df["lineage"].str.startswith( "B.1." ),"group"] = "B.1.X"

    df.loc[df["lineage"] == "B.1.1","group"] = "B.1.1"
    df.loc[df["lineage"].str.startswith(( "B.1.1.", "C", "P", "N", "R", "I", "L", "D", "M", "K", "J" )),"group"] = "B.1.1"

    df.loc[df["lineage"] == "B.1.2","group"] = "B.1.2"
    df.loc[df["lineage"].str.startswith( "B.1.2." ),"group"] = "B.1.2"

    df.loc[df["lineage"].isin(["B.1.429", "B.1.427"]), "group"] = "CA"

    return df


def split_clades( args ):
    lineages = load_clades( args.lineages )

    print( "Subsampled alignment is being split into the following clades." )
    print( lineages["group"].value_counts() )

    for clade, clade_md in lineages.groupby( "group" ):
        with tempfile.NamedTemporaryFile( suffix=".txt", mode="w+", delete=False ) as temp_clade:
            entries = clade_md["taxon"].to_list()
            entries.append( "CHN/Hubei-Wuhan/MN908947.3/2019-12-26" )
            entries.append( "CHN/Hubei-Wuhan/LR757998.1/2019-12-26" )
            entries = list( set( entries ) )
            temp_clade.write( "\n".join( entries ) )

        command = ["module load seqtk",
                   "seqtk subseq {} {} > {}".format( args.alignment, temp_clade.name, os.path.join( args.output, "clade_{}.fasta".format( clade ) ) )]

        call( " && ".join( command ), shell=True )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Splits alignment into clades based on pangolin lineage assignments" )

    # Initialize optional arguments
    parser.add_argument( "-a", "--alignment", help="location" )
    parser.add_argument( "-l", "--lineages", help="location" )
    parser.add_argument( "-o", "--output", help="location" )

    arguments = parser.parse_args()
    split_clades( arguments )
