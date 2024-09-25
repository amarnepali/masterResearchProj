#!/bin/bash
#SBATCH --ntasks=1 # indicate the number of core required
#SBATCH --job-name=nepa0034nsga3_dynamic200P
#SBATCH --output=/home/$FAN/%x-%j.out.txt
#SBATCH --error=/home/$FAN/%x-%j.err.txt
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --gres="gpu:1"
#SBATCH --time=5-0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nepa0034@flinders.edu.au

####################

source /home/nepa0034/.bashrc
conda activate PyEnv39

python3 optimization_dynamic_cm_200.py
