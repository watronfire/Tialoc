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


rule filter:
    message:
        """
        Filtering to
          - excluding strains in {params.exclude}
          - minimum genome length of {params.min_length}
          - Allowed date range of {params.min_date} -> {params.max_date}
        """
    group: "preparation"
    input:
        sequences = rules.combine_data.output.sequences,
        metadata = rules.combine_data.output.metadata
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
    group: "preparation"
    input:
        sequences = rules.filter.output.sequences,
    output:
        alignment = os.path.join( config["output"], "aligned.fasta" )
    params:       
        alignment_folder = os.path.join( config["output"], "alignment/" ),
        temp_alignment = os.path.join( config["output"], "alignment/filtered.fasta.aln" )
    threads: 16
    shell:
        """
        ViralMSA.py \
            -s {input.sequences} \
            -r {config[reference]} \
            -e {config[email]} \
            -o {params.alignment_folder} \
            -t {threads} &&
        mv {params.temp_alignment} {output.alignment}
        """


rule mask_vcf:
    message:
        """
        Mask bases in alignment based on provided VCF. Based on VCF provided by De Maio et al. at https://virological.org/t/masking-strategies-for-sars-cov-2-alignments/
        """
    group: "preparation"
    input:
        alignment = rules.align.output.alignment,
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
    group: "tree"
    output:
          collapsed_tree = os.path.join( config["output"], "tree/collapsed_tree.nwk" )
    shell:
         """
         {python} workflow/scripts/collapse_polytomies.py \
            --limit {config[collapse_polytomies][limit]} \
            --output {output.collapsed_tree} \
            {config[collapse_polytomies][tree_url]}
         """
