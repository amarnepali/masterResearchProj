#!/bin/bash
#SBATCH --ntasks=1 # indicate the number of core required
#SBATCH --job-name=nepa0034orbOpt
#SBATCH --output=orbOpt-log.txt
#SBATCH --error=logs/%A.error
#SBATCH --partition=melfu
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --gres="gpu:0"
#SBATCH --time=1-0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nepa0034@flinders.edu.au

####################

source /home/nepa0034/.bashrc
conda activate PyEnv39

python3 optimization_nsga3_simple.py
