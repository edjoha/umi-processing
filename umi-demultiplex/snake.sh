#!/bin/bash

#SBATCH --job-name=demux
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --output=/fast/users/altwassr_c/scratch/slurm_logs/demux-%x.%j.out
#SBATCH --error=/fast/users/altwassr_c/scratch/slurm_logs/demux-%x.%j.err

echo $(date -u) 'Start'
snakemake \
    --nt \
    --jobs 40 \
    --keep-going \
    --rerun-incomplete \
    --workflow-profile=my_profile \
    -p --use-conda --conda-frontend mamba --conda-prefix /data/cephfs-1/work/projects/damm-targetseq/conda \
    --cores 40
echo $(date -u) 'Finished' #time?

#my_profile is basically the same as cubi-v1, but with increased number of restarts and latency-wait time


    #--dry-run \
    # --restart-times 2 \
    # --reason \
