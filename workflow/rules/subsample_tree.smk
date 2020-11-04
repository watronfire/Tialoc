import pandas as pd
from subprocess import call

rule parse_llama_output:
    message: "subsamples sequences from llama output catchment files"
    group: "subsample_tree"
    input:
        llama_output = directory( os.path.join( config["output"], "llama_output" ) ),
        metadata = rules.add_interest.output.metadata
    output:
        sequences = os.path.join( config["output"], "selected_sequences.txt" )
    shell:
        """
        {python} workflow/scripts/parse_llama_output.py \
            --metadata {input.metadata} \
            --llama {input.llama_output} \
            --output {output.sequences}
        """

rule extract_llama_output:
    message: "Filter metadata and alignment to sequences specified by llama output."
    group: "subsample_tree"
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

        md = pd.read_csv( input.metadata )
        md = md.loc[md["strain"].isin( seqs )]
        print( "{} of {} entries found in metadata".format( len( md ), len( seqs ) ) )
        md.to_csv( output.subsampled_metadata, index=False )

        # Filter alignment
        command = ["module load seqtk",
                   "seqtk subseq {} {} > {}".format( input.alignment, input.sequences, output.subsampled_alignment )]
        call( " && ".join( command ), shell=True )