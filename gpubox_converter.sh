#!/bin/bash -l
#SBATCH --job-name="gpubox_job"
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --output=gpubox_job.%j.o
#SBATCH --error=gpubox_job.%j.e
#SBATCH --mem=120G
#SBATCH --export=ALL

export TMPDIR=$home2/tmp
export TMP=$home2/tmp
export TEMP=$home2/tmp


python /home/nbarry/MWA/MWA_data_analysis/gpubox_converter.py
