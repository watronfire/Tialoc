#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q gpu -o /gpfs/home/natem/logs/tialoc.txt -j oe

conda deactivate
module load python/3.6.3
module load R

export AUGUR_RECURSION_LIMIT=3000

cd /gpfs/home/natem/scripts/Tialoc
snakemake -k -j 50 \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --cluster-config config/cluster.json \
  --cluster "qsub -V -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
