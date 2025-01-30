#!/bin/bash
#SBATCH --job-name=Pstar
#SBATCH --partition=gpp
#SBATCH --time=11:59:00
#SBATCH --output=output/out
#SBATCH --error=output/err
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mail-type=all
#SBATCH --mail-user=h.sanchez@ua.pt

srun ./bhac -i amrvac_start.par

