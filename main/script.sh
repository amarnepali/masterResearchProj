#!/bin/bash

SBATCH --job-name=orbOpt
SBATCH --output=orbOpt-log.txt
SBATCH --ntasks=1
SBATCH --cpus-per-task=32
SBATCH --mem=32GB
SBATCH --time=12:00:00

####################

# conda init bash
conda activate PyEnv39

python3 ConstellationProj/masterResearchProj/main/Optimization.py