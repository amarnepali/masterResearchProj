#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24 # indicate the number of core required
#SBATCH --partition=amilan
#SBATCH --job-name=nepa0034orbOpt
#SBATCH --output=orbOpt-log.txt
#SBATCH --error=logs/%A.error
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nepa0034@flinders.edu.au

####################

source /home/nepa0034/.bashrc
conda activate PyEnv39

python3 Optimization.py
