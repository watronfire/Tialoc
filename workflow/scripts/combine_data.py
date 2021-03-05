import pandas as pd
from Bio import SeqIO
import os
import argparse
import unicodedata
import re

used_names = list()
count = 1

def get_SEARCH_name( entry, dedup=False, remove_incomplete=True ):
    sname = list()
    sname.append( entry["country_abv"] )
    locations = entry[["division","location"]].to_list()
    locations = [str(j) for j in locations]
    if "nan" in locations: locations.remove( "nan" )
    if "?" in locations: locations.remove( "?" )
    sname.append( "-".join( locations ) )
    sname.append( entry["accession"] )
    sname.append( str( entry["date"] ) )
    new_name = "/".join( sname )
    new_name = new_name.replace( " ", "" )
    new_name = new_name.replace( "'", "" )
    if dedup:
        if new_name in used_names:
            print( "{} already used".format( new_name ) )
            return "bad_name"
        used_names.append( new_name )

    if remove_incomplete:
        if entry["country_abv"] == "?":
            return "bad_name"

    return new_name

def load_county_dict( file_location ):
    countries = pd.read_csv( file_location )
    countries = countries.set_index( "Name" )
    return countries["Alpha-3"].to_dict()

def parse_exclude( file_location ):
    return_list = list()
    with open( file_location, "r" ) as exluded:
        for line in exluded:
            if not line.startswith( "#" ) and line != "\n":
                return_list.append( line.strip() )
    return return_list

def identify_gisaid_files( loc ):
    md_regex = re.compile( "metadata_[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}_[0-9]{1,2}-[0-9]{1,2}\.tsv" )
    seq_regex = re.compile( "sequences_[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}_[0-9]{1,2}-[0-9]{1,2}\.fasta" )
    for j in os.listdir( loc ):
        if md_regex.match( j ):
            metadata = os.path.join( loc, j )
        elif seq_regex.match( j ):
            sequences = os.path.join( loc, j )
    assert 'metadata' in locals(), "Unable to find GISAID metadata. Either file is not in the directory or regex is incorrect."
    assert 'sequences' in locals(), "Unable to find GISAID sequences. Either file is not in the directory or regex is incorrect."

    return sequences, metadata

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser( 
            description="Combine metadata and sequences from SEARCH repository and GISAID",
            formatter_class= argparse.ArgumentDefaultsHelpFormatter
    )

    # TODO gseqs and gmetadata arguments aren't necessary if we place them in the same location every time.
    parser.add_argument( "--search", type=str, required=True, help="SEARCH repository" )
    parser.add_argument( "--country", type=str, required=True, help="Path to country dictionary in csv format" )
    parser.add_argument( "--exclude", type=str, required=True, help="Path to excluded sequence list" )

    args = parser.parse_args()

    cd = load_county_dict( args.country )

    excluded_list = parse_exclude( args.exclude )
    gseqs, gmetadata = identify_gisaid_files( args.search )

    print(  "Formating gisaid metadata... ", end="" )
    gisaid_md = pd.read_csv( gmetadata, delimiter="\t" )
    gisaid_md = gisaid_md.loc[~gisaid_md["strain"].isin( excluded_list)]
    gisaid_md["accession"] = gisaid_md["genbank_accession"]
    gisaid_md.loc[gisaid_md["genbank_accession"]=="?","accession"] = gisaid_md["strain"].apply( lambda x : x.replace( "/", "" ) )
    gisaid_md["country_abv"] = gisaid_md["country"].map( cd )

    gisaid_md["name"] = gisaid_md["strain"]
    gisaid_md = gisaid_md.fillna( "?" )

    # Remove non-ascii characters where present.
    for i in gisaid_md.columns:
        try:
            gisaid_md[i] = gisaid_md[i].apply( lambda val: unicodedata.normalize( 'NFKD', val ).encode( 'ascii', 'ignore' ).decode() )
        except TypeError:
            continue

    gisaid_md["location_str"] = "-" + gisaid_md["location"].astype(str)
    gisaid_md.loc[gisaid_md["location_str"].str.contains("\?"),"location_str"] = ""
    gisaid_md["location_str"] = gisaid_md["division"].astype(str) + gisaid_md["location_str"].astype(str)
    gisaid_md["strain"] = gisaid_md["country_abv"] \
                          + "/" + gisaid_md["location_str"] \
                          + "/" + gisaid_md["accession"] \
                          + "/" + gisaid_md["date"]
    gisaid_md["strain"] = gisaid_md["strain"].str.replace( " ", "" )
    gisaid_md["strain"] = gisaid_md["strain"].str.replace( "'", "" )

    gisaid_md["strain"] = gisaid_md["strain"].apply(lambda val: unicodedata.normalize('NFKD', val).encode('ascii', 'ignore').decode())

    # Reorder according to nextstrain order
    gisaid_md = gisaid_md[["strain",
                           "virus",
                           "gisaid_epi_isl",
                           "genbank_accession",
                           "date",
                           "region",
                           "country",
                           "division",
                           "location",
                           "region_exposure",
                           "country_exposure",
                           "division_exposure",
                           "segment",
                           "length",
                           "host",
                           "age",
                           "sex",
                           "originating_lab",
                           "submitting_lab",
                           "authors",
                           "url",
                           "title",
                           "date_submitted",
                           "name",
                           "pangolin_lineage"]]

    gisaid_md = gisaid_md.set_index( "name" )
    
    with open( os.path.join( args.search, "metadata.csv" ), "r", encoding="ascii", errors="ignore" ) as temp_md:
        search_md = pd.read_csv( temp_md )

    # Determine Country, division, and location
    search_md[["country","division", "location"]] = search_md["location"].str.split( "/", expand=True )
    search_md.loc[search_md["country"]=="MEX","country"] = "Mexico"

    # cleanup location
    search_md.loc[search_md["location"].isnull(),"location"] = "?"

    # cleanup country
    search_md["country_abv"] = search_md["country"].replace( cd )

    # Determine Region
    search_md["region"] = "North America"
    search_md.loc[search_md["country"]=="Jordan","region"] = "Asia"

    search_md = search_md.rename(columns={"ID":"accession", 
                                          "collection_date": "date", 
                                          "gb_accession" : "genbank_accession", 
                                          "gisaid_accession" : "gisaid_epi_isl" } )
    search_md["location_str"] = "-" + search_md["location"].astype(str)
    search_md.loc[search_md["location_str"].str.contains("\?"),"location_str"] = ""
    search_md["location_str"] = search_md["division"].astype(str) + search_md["location_str"].astype(str)
    search_md["strain"] = search_md["country_abv"] \
                          + "/" + search_md["location_str"] \
                          + "/" + search_md["accession"] \
                          + "/" + search_md["date"]
    search_md["strain"] = search_md["strain"].str.replace( " ", "" )
    search_md["strain"] = search_md["strain"].str.replace( "'", "" )
    search_md["ID"] = search_md["accession"]
    search_md = search_md.set_index( "accession" ) 

    # Add missing data
    search_md["virus"] = "ncov"
    search_md["region_exposure"] = search_md["region"]
    search_md["country_exposure"] = search_md["country"]
    search_md["division_exposure"] = search_md["division"] 
    search_md["segment"] = "genome"
    search_md["host"] = "Human"
    search_md["age"] = "?"
    search_md["sex"] = "?"
    search_md["submitting_lab"] = "Andersen lab at Scripps Research"
    search_md["url"] = "https://github.com/andersen-lab/HCoV-19-Genomics"
    search_md["title"] = "?"
    search_md["date_submitted"] = search_md["date"]
    search_md["pangolin_lineage"] = "?"

    # We'll initialize this later
    search_md["length"] = 0

    #Rename columns and reorder according to GISAID 
    search_md = search_md[["strain",
                           "virus",
                           "gisaid_epi_isl",
                           "genbank_accession",
                           "date",
                           "region",
                           "country",
                           "division",
                           "location",
                           "region_exposure",
                           "country_exposure",
                           "division_exposure",
                           "segment",
                           "length",
                           "host",
                           "age",
                           "sex",
                           "originating_lab",
                           "submitting_lab",
                           "authors",
                           "url",
                           "title",
                           "date_submitted",
                           "pangolin_lineage"]]

    search_md = search_md.fillna( "?" )

    seq_db = search_md["strain"].to_dict()
    gisaid_db = search_md["gisaid_epi_isl"].to_dict()
    search_location = os.path.join( args.search, "consensus_sequences" )
    files = [i for i in os.listdir(search_location) if i.lower().endswith((".fa",".fasta"))]

    search_seqs = list()
    search_ids = list()
    len_dict = dict()
    for file in files:
        name = file.split( "." )[0]
        seq = SeqIO.read( os.path.join( search_location, file ), "fasta" )
        seq.id = seq_db[name]
        search_ids.append( gisaid_db[name] )
        seq.description = ""
        seq.name = ""
        search_seqs.append( seq )
        len_dict[name] = len( seq.seq )
    search_md["length"] = search_md.index.map( len_dict )

    gisaid_seqs = SeqIO.parse( gseqs, "fasta" )

    count = 0
    exclude_count = 0
    strain_db = gisaid_md["strain"].to_dict()
    gisaid_db = gisaid_md["gisaid_epi_isl"].to_dict()
    pangolin_db = gisaid_md["pangolin_lineage"].to_dict()
    pangolin_dict = dict()
    for seq in gisaid_seqs:
        if seq.name in excluded_list:
            exclude_count += 1
            continue        

        if seq.name not in strain_db.keys():
            print( "Unable to find: {}" )
            continue

        if gisaid_db[seq.name] in search_ids:
            pangolin_dict[gisaid_db[seq.name]] = pangolin_db[seq.name]
            count += 1
            continue

        seq.id = strain_db[seq.name]
        seq.description = ""
        seq.name = ""
        search_seqs.append( seq )
    print( "Done" )

    print( "\t{} sequences already in dataset.".format( count ) )
    print( "\t{} sequences removed due to exclude list.".format( exclude_count ) )

    # Add pangolin lineage from gisaid if available
    search_md["pangolin_lineage"] = search_md["gisaid_epi_isl"].map( pangolin_dict )

    gisaid_md = gisaid_md.loc[gisaid_md["strain"].isin([i.id for i in search_seqs])]

    combined_md = pd.concat( [search_md.reset_index( drop=True ), gisaid_md], ignore_index=True )

    # This could be wrong
    combined_md = combined_md.groupby( "strain" ).first().reset_index()

    # Add search Metadata
    combined_md["search"] = "Other"
    combined_md.loc[combined_md["country"] == "Mexico","search"] = "Mexico"
    combined_md.loc[combined_md["country"] == "USA", "search"] = "USA"
    combined_md.loc[combined_md["division"] == "BajaCalifornia", "search"] = "Baja California"
    combined_md.loc[combined_md["division"] == "Baja California", "search"] = "Baja California"
    combined_md.loc[combined_md["division"] == "California", "search"] = "California"
    combined_md.loc[combined_md["location"].str.startswith( "San Diego" ), "search"] = "San Diego"
    combined_md.loc[combined_md["division"] == "Texas", "search"] = "Texas"
    combined_md.loc[combined_md["division"] == "Arizona", "search"] = "Arizona"
    combined_md.loc[combined_md["division"] == "New Mexico", "search"] = "New Mexico"

    combined_md.to_csv( os.path.join( args.search, "metadata.tsv" ), index=False, sep="\t", encoding="ascii" )

    SeqIO.write( search_seqs, os.path.join( args.search, "sequences.fasta" ), "fasta" )

