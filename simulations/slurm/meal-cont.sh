#!/bin/bash
# Job name:
#SBATCH --job-name=cont-sim
#
# Processors (1 node = X cores):
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=X
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=your.email@here.com
#
# Job output:
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#
## Command(s) to run:
cd ~/pub_temvip-framework/simulations
#
## Run the script
R CMD BATCH --no-save --no-restore \
  R/meal_obs-continuous.R \
  logs/meal_obs-continuous.Rout
