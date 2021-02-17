import datetime

rule combine_data:
    message: "Combine data from GISAID and SEARCH github repository, while also renaming sequences to match useful format"
    group: "preparation"
    input:
        gisaid_seqs = config["gisaid_seqs"],
        gisaid_metadata = config["gisaid_md"],
        search_md = rules.download_search_data.output.metadata,
        exclude = rules.download_exclude.output.exclude
    params:
        search_repo = os.path.join( config["output"], "HCoV-19-Genomics/" ),
        country_dict = config["combine_data"]["country_codes"],
    output:
        sequences = os.path.join( config["output"], "HCoV-19-Genomics/sequences.fasta" ),
        metadata = os.path.join( config["output"], "HCoV-19-Genomics/metadata.tsv" )
    shell:
         """
         {python} workflow/scripts/combine_data.py \
            --search {params.search_repo} \
            --gseqs {input.gisaid_seqs} \
            --gmetadata {input.gisaid_metadata} \
            --country {params.country_dict} \
            --exclude {input.exclude}
         """


rule add_interest:
    message: "Add focus to metadata"
    group: "preparation"
    input:
        metadata = rules.combine_data.output.metadata
    output:
        metadata = os.path.join( config["output"], "HCoV-19-Genomics/metadata_interest.tsv" )
    shell:
        """
        {python} workflow/scripts/add_interest.py \
            --metadata {input.metadata} \
            --interest {config[focus]:q} \
            --output {output.metadata}
        """


rule filter:
    message:
        """
        Filtering to
          - excluding strains in {params.exclude}
          - minimum genome length of {params.min_length}
          - Allowed date range of {params.min_date} -> {params.max_date}
        """
    input:
        sequences = rules.combine_data.output.sequences,
        metadata = rules.add_interest.output.metadata
    params:
        include = config["include"],
        exclude = config["exclude"],
        min_length = config["filter"]["min_length"],
        exclude_where = config["filter"]["exclude_where"],
        group_by = config["filter"]["group_by"],
        sequences_per_group = config["filter"]["sequences_per_group"],
        min_date = config["filter"]["min_date"],
        max_date = datetime.datetime.today().strftime( "%Y-%m-%d" )
    output:
        sequences = os.path.join( config["output"], "filtered.fasta" )
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {params.include} \
            --exclude {params.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output.sequences}
        """


rule align:
    message:
        """
        Aligning sequences to {config[reference]}
          - gaps relative to reference are considered real
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = config["reference"]
    output:
        sam_file = temp( os.path.join( config["output"], "aligned.sam" ) )
    threads: 16
    shell:
        """
        minimap2 \
            -a -x asm5 \
            -t {threads} \
            {input.reference} \
            {input.sequences} \
            -o {output.sam_file}
        """

rule sam_2_fasta:
    input:
        sam = rules.align.output.sam_file,
        reference = config["reference"]
    message: "Running datafunk to trim and pad against the reference"
    params:
        trim_start = config["align"]["trim_start"],
        trim_end = config["align"]["trim_end"]
    output:
        alignment = os.path.join(config["output"],"aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.alignment} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts 
        """


rule mask:
    message:
        """
        Mask bases in alignment based on provided VCF. Based on VCF provided by De Maio et al. at https://virological.org/t/masking-strategies-for-sars-cov-2-alignments/
        """
    input:
        alignment = rules.sam_2_fasta.output.alignment,
        mask = rules.download_mask.output.vcf
    output:
        alignment = os.path.join( config["output"], "masked.fasta" )
    shell:
        """
        {python:q} workflow/scripts/mask_alignment_vcf.py \
            -i {input.alignment} \
            -o {output.alignment} \
            -v {input.mask}
        """


rule collapse_polytomies:
    message: "Collapse edges of tree that are indistinguishable from polytomies"
    group: "tree"
    output:
        collapsed_tree = os.path.join( config["output"], "tree/collapsed_tree.nwk" )
    params:
        # --url {config[collapse_polytomies][tree_url]}
        tree = os.path.join( config["output"], "tree_whole.newick" )
    shell:
         """
         {python} workflow/scripts/collapse_polytomies.py \
            --limit {config[collapse_polytomies][limit]} \
            --output {output.collapsed_tree} \
            --path {params.tree}
            --rename {config[gisaid_md]}
         """

rule prune_tree:
    message: "Prunes input tree to match alignment, metadata, and focus"
    group: "align"
    input:
        tree = rules.collapse_polytomies.output.collapsed_tree,
        alignment = rules.mask.output.alignment,
        metadata = rules.add_interest.output.metadata
    output:
        global_tree = os.path.join( config["output"], "data-dir/global.tree" ),
        global_alignment = os.path.join( config["output"], "data-dir/alignment.fasta" ),
        global_metadata = os.path.join( config["output"], "data-dir/metadata.csv" ),
        query_alignment = os.path.join( config["output"], "data-dir/query.fasta" ),
        query_metadata = os.path.join( config["output"], "data-dir/query.csv" )
    params:
        outdir = os.path.join( config["output"], "data-dir" )
    shell:
        """
        {python} workflow/scripts/prune_tree.py \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --outdir {params.outdir}
        """


rule generate_llama_script:
    message: "Generating cluster script for llama."
    group: "align"
    input:
        global_tree = rules.prune_tree.output.global_tree,
        global_alignment = rules.prune_tree.output.global_alignment,
        global_metadata = rules.prune_tree.output.global_metadata,
        query_alignment = rules.prune_tree.output.query_alignment,
        query_metadata = rules.prune_tree.output.query_metadata
    output:
        llama_script = os.path.join( config["output"], "llama_script.sh" )
    shell:
        "{python} workflow/scripts/generate_llama_script.py --outdir {config[output]}"
