rule download_search_data:
    message: "Download data from SEARCH github repository"
    group: "init"
    output:
        metadata = os.path.join( config["output"], "HCoV-19-Genomics/metadata.csv" )
    shell:
        """
        cd {config[output]} 
        git clone https://github.com/andersen-lab/HCoV-19-Genomics.git
        """

rule download_exclude:
    message: "Updating excluded sequence list from ncov repository"
    group: "init"
    output:
        exclude = config["exclude"]
    shell:
        """
        curl https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt \
            --output {output.exclude}
        """

rule download_mask:
    message: "Collect vcf mask from repository"
    output:
        vcf = "resources/mask.vcf"
    shell:
        """
        curl https://raw.githubusercontent.com/watronfire/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
            --output {output.vcf}
        """