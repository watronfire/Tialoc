#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --partition=shared
#SBATCH --output=/gpfs/home/natem/logs/tialoc.txt

module load python/3.8.3
module load R

export AUGUR_RECURSION_LIMIT=3000

cd $SLURM_SUBMIT_DIR
echo $pwd
snakemake -k -j 50 \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --cluster-config config/cluster.json \
  --cluster "sbatch --time={cluster.walltime} --mem={cluster.mem} -c {cluster.n} --partition={cluster.queue} --output={cluster.logfile}" /gpfs/home/natem/analysis/2021.02.16_hcov/tree/collapsed_tree.nwk
