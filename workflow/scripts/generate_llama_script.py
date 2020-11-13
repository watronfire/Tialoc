import argparse
import os

script_template = """
#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=16 -l mem=64gb -q workq -o /gpfs/home/natem/logs/llama.txt -j oe

module rm python/3.6.3
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

source ~/miniconda3/etc/profile.d/conda.sh
conda activate llama

cd {outdir}
llama \
    --input data-dir/query.csv \
    --datadir data-dir/ \
    --input-column strain \
    --data-column strain \
    --node-summary country \
    --no-temp \
    -nr \
    --fasta data-dir/query.fasta \
    --threads 16 \
    --outdir llama_output/
"""

def generate_llama_script( args ):
    with open( os.path.join( args.outdir, "llama_script.sh" ), "w" ) as script:
        script.write( script_template.format( outdir=args.outdir ) )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Generates a cluster script for llama as it can't be executed by snakemake at the moment." )
    parser.add_argument( "-o", "--outdir", help="output directory" )

    arguments = parser.parse_args()
    
    generate_llama_script( arguments )
