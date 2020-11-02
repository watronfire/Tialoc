
rule run_llama:
    message: "Pull substrees from tree corresponding to sequences in query metadata"
    group: "subsample"
    input:
        global_tree = rules.prune_tree.output.global_tree,
        global_alignment = rules.prune_tree.output.global_alignment,
        global_metadata = rules.prune_tree.output.global_metadata,
        query_alignment = rules.prune_tree.output.query_alignment,
        query_metadata = rules.prune_tree.output.query_metadata
    output:
        collapsed_nodes = os.path.join( config["output"], "llama_output/catchment_trees/localcollapsed_nodes.csv" ),
        local_trees = directory( os.path.join( config["output"], "llama_output/local_trees" ) )
    params:
        data_dir = os.path.join( config["output"], "data-dir" ),
        output_dir = os.path.join( config["output"], "llama_output" )
    conda: "/gpfs/home/natem/scripts/Tialoc/workflow/envs/llama.yaml"
    threads: 16
    shell:
        """
        llama \ 
            --input {input.query_metadata} \
            --datadir {params.data_dir} \
            --input-column {config[run_llama][input_column]} \
            --data-column {config[run_llama][input_column]} \
            --node-summary {config[run_llama][summarize]} \
            --fasta {input.query_alignment} \
            --outdir {params.output_dir} \
            --threads {threads}
            -nr
        """
