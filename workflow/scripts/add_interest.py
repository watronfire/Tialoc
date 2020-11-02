import argparse
import pandas as pd

def parse_interest( interest_str ):
    return_list = list()
    for i in interest_str.split( ";" ):
        resolution, values = i.split( ":", 1 )
        for j in values.split( "," ):
            return_list.append( (resolution, j) )
    return return_list

def main( args ):
    rejected_values = ["Cruise_Ship_2", "Cruise_Ship_1", "Grand Princess Cruise Ship"]

    interest = parse_interest( args.interest )

    df = pd.read_csv( args.metadata, sep="\t" )

    df["interest"] = "non-interest"
    for resolution, value in interest:
        df.loc[df[resolution] == value, "interest"] = "interest"

    df.loc[df["location"].isin( rejected_values ), "interest"] = "non-interest"
    df["lineage"] = "?"

    df.to_csv( args.output, sep="\t", index=False )
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Add focus to metadata" )

    # Initialize optional arguments
    parser.add_argument( "-o", "--output", help="ouput location" )
    parser.add_argument( "-i", "--interest", help="Focus to add. In the format resolution1:value1,value2;resolution2:value3,value4" )
    parser.add_argument( "-m", "--metadata", help="Metadata to modify" )

    arguments = parser.parse_args()

    main( arguments )
