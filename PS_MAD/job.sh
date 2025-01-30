#!/bin/bash
#SBATCH --job-name=PstarMAD
#SBATCH --partition=gpp
#SBATCH --time=71:59:00
#SBATCH --output=output/out
#SBATCH --error=output/err
#SBATCH --ntasks=224
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mail-type=all
#SBATCH --mail-user=jose.olvera-meneses@uni-tuebingen.de


srun ./bhac -restart 618 -slice 619 -collapse 619

