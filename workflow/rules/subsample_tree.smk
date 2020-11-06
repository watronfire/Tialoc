import pandas as pd
from subprocess import call

rule parse_llama_output:
    message: "subsamples sequences from llama output catchment files"
    input:
        closes_in_db = os.path.join( config["output"], "llama_output/closest_in_db.csv" ),
        metadata = expand( os.path.join( config["output"], "data-dir/{md}.csv" ), md=["metadata", "query"] ) 
    params:
        llama_output = directory( os.path.join( config["output"], "llama_output" ) )
    output:
        sequences = os.path.join( config["output"], "selected_sequences.txt" )
    shell:
        """
        {python} workflow/scripts/parse_llama_output.py \
            --metadata {input.metadata} \
            --llama {params.llama_output} \
            --output {output.sequences}
        """

rule extract_llama_output:
    message: "Filter metadata and alignment to sequences specified by llama output."
    input:
        alignment = rules.mask.output.alignment,
        metadata = rules.add_interest.output,
        sequences = rules.parse_llama_output.output.sequences
    output:
        subsampled_alignment = os.path.join( config["output"], "output/subsampled_alignment.fasta" ),
        subsampled_metadata = os.path.join( config["output"], "output/subsampled_metadata.csv" )
    run:
        # load sequences
        with open( input.sequences, "r" ) as seq_file:
            seqs = [line.strip() for line in seq_file]
        seqs.extend( config["extract_llama_output"]["include"].split( " " ) )

        md = pd.read_csv( str( input.metadata ), sep="\t" )
        md = md.loc[md["strain"].isin( seqs )]
        print( "{} of {} entries found in metadata".format( len( md ), len( seqs ) ) )
        md.to_csv( output.subsampled_metadata, index=False )

        # Filter alignment
        command = ["module load seqtk",
                   "seqtk subseq {} {} > {}".format( input.alignment, input.sequences, output.subsampled_alignment )]
        call( " && ".join( command ), shell=True )


rule build_tree:
    message: "Generate tree from llama subsampling"
    input:
        alignment = rules.extract_llama_output.output.subsampled_alignment
    output:
        tree = os.path.join( config["output"], "output/subsampled_tree.newick" ),
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --method "iqtree" \
            --output {output.tree} \
            --substitution-model {config[build_tree][model]} \
            --nthreads 16
        """

