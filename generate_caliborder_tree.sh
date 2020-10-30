#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q gpu -o /gpfs/home/natem/logs/nextstrain.txt -j oe

#conda deactivate
#module load python/3.6.3

#cd /Users/natem/Dropbox (Scripps Research)/Personal/Code/Python/caliborder_generator/workflow
snakemake -k -j 50 \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml