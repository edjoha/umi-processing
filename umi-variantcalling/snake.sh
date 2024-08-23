#!/bin/bash

# export SNAKEMAKE_SLURM_DEBUG=1

#SBATCH --job-name=variantcalling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=24000M
#SBATCH --output=/data/cephfs-1/work/projects/damm-belove/scratch/A4314_reseq/slurm/%x.%j.out
#SBATCH --error=/data/cephfs-1/work/projects/damm-belove/scratch/A4314_reseq/slurm/%x.%j.err

snakemake \
    --jobs 10 \
    --cluster-config /data/gpfs-1/users/cofu10_c/work/pipelines/umi-processing/config/cluster_config.yaml \
    --profile=cubi-v1 \
    --restart-times 2 \
    --keep-going \
    --rerun-incomplete \
    --use-conda \
    --conda-prefix=/data/gpfs-1/users/cofu10_c/scratch/A4314_reseq/envs     
    # --verbose \
    #--touch \
    #--skip-script-cleanup \
    #--reason 
    #--until annovar
    #--until table_to_anno \
