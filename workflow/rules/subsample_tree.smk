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

        md = pd.read_csv( str( input.metadata ), sep="\t" )
        md = md.loc[md["strain"].isin( seqs )]
        print( "{} of {} entries found in metadata".format( len( md ), len( seqs ) ) )
        md.to_csv( output.subsampled_metadata, index=False )

        # Filter alignment
        command = ["module load seqtk",
                   "seqtk subseq {} {} > {}".format( input.alignment, input.sequences, output.subsampled_alignment )]
        call( " && ".join( command ), shell=True )


rule identify_lineages:
    message: "Use pangolin to identify lineages"
    input:
        alignment = rules.extract_llama_output.output.subsampled_alignment
    params:
        outdir = os.path.join( config["output"], "output/" )
    conda:
        "../envs/pangolin.yaml"
    output:
        lineages = os.path.join( config["output"], "output/lineage_report.csv" )
    shell:
        """
        pangolin {input.alignment} -o {params.outdir} -t 8
        """


rule split_lineages:
    input:
        alignment = rules.extract_llama_output.output.subsampled_alignment,
        lineages = rules.identify_lineages.output.lineages
    params:
        outdir = os.path.join( config["output"], "clade_alignment" )
    output:
        clade_alignments = expand( os.path.join( config["output"], "clade_alignment/clade_{clade}.fasta" ), clade=config["clades"] )
    shell:
        """
        {python} workflow/scripts/split_clade_alignments.py \
            --alignment {input.alignment} \
            --lineages {input.lineages} \
            --clades {config[clades]} \
            --output {params.outdir}
        """


rule build_clade_tree_fast:
    message: "Generate treef rom llama subsampling: {wildcards.clade}"
    input:
        alignment = os.path.join( config["output"], "clade_alignment/clade_{clade}.fasta" )
    params:
        clade = "{wildcards.clade}"
    output:
        tree = os.path.join( config["output"], "clade_trees/subsampled_{clade}_tree.unrooted.newick" )
    threads: 16
    log:
        os.path.join( config["output"], "logs/build_tree_{clade}.log" )
    shell:
        """
        echo "{input.alignment} {params.clade}" &&
        export OMP_NUM_THREADS={threads} &&
        module load fasttree &&
        fasttreemp -nosupport -nt {input.alignment} > {output.tree} 2> {log}
        """

rule root_tree:
    input:
        tree = rules.build_clade_tree_fast.output.tree,
    params:
        clade = "{wildcard.clade}",
        outgroup = config["build_tree"]["outgroup"]
    output:
        tree = os.path.join( config["output"], "clade_trees/subsampled_{clade}_tree.newick" ),
    log:
        os.path.join( config["output"], "logs/root_tree_{clade}.log" ),
    shell:
        """
        echo "{params.clade} {params.outgroup}"
        jclusterfunk reroot \
          --format newick \
          -i {input.tree} \
          -o {output.tree} \
          --outgroup {params.outgroup} &>> {log}
        """

rule build_whole_tree:
    message: "Generate tree from llama subsampling"
    input:
        alignment = rules.extract_llama_output.output.subsampled_alignment
    params:
        outgroup = config["build_tree"]["outgroup"].replace( "/", "_X_X_" )
    output:
        tree = os.path.join( config["output"], "output/ml_trees/subsampled_tree.newick" )
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --method "iqtree" \
            --output {output.tree} \
            --substitution-model {config[build_tree][model]} \
            --nthreads 16 \
            --tree-builder-args='-o {params.outgroup}'
        """

rule collapse_polytomies_alt:
    message: "Collapse polytomies in IQTree output: {wildcards.clade}."
    input:
        tree = rules.root_tree.output.tree
    output:
        collapsed_tree = os.path.join( config["output"], "clade_trees/subsampled_{clade}_tree_collapsed.newick" )
    shell:
         """
         {python} workflow/scripts/collapse_polytomies.py \
            --limit {config[collapse_polytomies][limit]} \
            --output {output.collapsed_tree} \
            --path {input.tree}
         """

rule clock_rate_filter:
    message: "remove tips from tree which violate infered or specified clock rate: {wildcards.clade}"
    input:
        tree = rules.collapse_polytomies_alt.output.collapsed_tree
    output:
        filtered_tree = os.path.join( config["output"], "clade_trees/subsampled_{clade}_tree_filtered.newick" )
    shell:
        """
        {python} workflow/scripts/clock_filter.py \
            --input {input.tree} \
            --clock-rate {config[clock_rate_filter][clock_rate]} \
            --iqd {config[clock_rate_filter][iqd]} \
            --output {output.filtered_tree}
        """


rule second_pruning:
    message: "Subject tips which were added after llama output parsing to pruning: {wildcards.clade}"
    input:
        tree = rules.clock_rate_filter.output.filtered_tree,
        metadata = rules.extract_llama_output.output.subsampled_metadata,
    output:
        pruned_tree = os.path.join( config["output"], "output/ml_trees/subsampled_{clade}_tree.newick" )
    shell:
        """
        {python} workflow/scripts/prune_interest.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.pruned_tree}
        """


rule resolve_polytomies:
    message: "Resolve output tree into a strickly bifurcating tree for beast: {wildcards.clade}"
    input:
        tree = rules.second_pruning.output.pruned_tree
    output:
        bi_tree = os.path.join( config["output"], "temp/subsampled_{clade}_tree.bifurcating.newick" )
    shell:
        """
        module load R &&
        Rscript workflow/scripts/resolve_polytomies.R {input.tree} {output.bi_tree} 
        """


rule generate_nexus:
    message: "Combine output alignment and bifuricating tree into nexus format for beast: {wildcards.clade}"
    input:
        tree = rules.second_pruning.output.pruned_tree,
        alignment = os.path.join( config["output"], "clade_alignment/clade_{clade}.fasta" ),
        metadata = rules.extract_llama_output.output.subsampled_metadata
    output:
        nexus_file = os.path.join( config["output"], "output/beast_input/clade_{clade}.nexus" ),
        traits = os.path.join( config["output"], "output/beast_input/clade_{clade}_traits.tsv" )
    shell:
        """
        {python} workflow/scripts/generate_nexus.py \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --fields search \
            --output {output.nexus_file} \

        """
