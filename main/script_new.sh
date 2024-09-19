#!/bin/bash
#SBATCH --job-name long_job_nsga3
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nepa0034@flinders.edu.au


source /home/nepa0034/.bashrc
conda activate PyEnv39

python3 optimization_nsga3_simple.py
