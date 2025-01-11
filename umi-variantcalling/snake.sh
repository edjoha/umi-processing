#!/bin/bash

# export SNAKEMAKE_SLURM_DEBUG=1

#SBATCH --job-name=variantcalling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=24000M
#SBATCH --output=/fast/users/altwassr_c/scratch/slurm_logs/%x.%j.out
#SBATCH --error=/fast/users/altwassr_c/scratch/slurm_logs/%x.%j.err

echo $(date -u) 'Start variant calling'
snakemake \
    --jobs 60 \
    --keep-going \
    --rerun-incomplete \
    --workflow-profile=my_profile \
    --cluster-config ~/work/umi-data-processing/config/cluster_config.yaml \
    -p --use-conda --conda-frontend mamba --conda-prefix /data/cephfs-1/work/projects/damm-targetseq/conda \
    --cores
    #-r \ ???
    #--nt \ >>> temp() ignored
echo $(date -u) 'Finished variant calling'

#write a log file??

# --touch \
# --skip-script-cleanup \
# --reason \

# --until annovar \
