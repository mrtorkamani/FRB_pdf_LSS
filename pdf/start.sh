#!/bin/bash
#SBATCH --partition=intelsr_short
#SBATCH --time=00:58:00
#SBATCH --output=secondslurm-%j.out
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1

python3 pipeline.py 
