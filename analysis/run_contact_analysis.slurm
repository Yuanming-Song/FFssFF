#!/bin/bash

#SBATCH --job-name=AA_contact_analysis
#SBATCH --account=dtobias_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --out=contact_analysis.out
#SBATCH --err=contact_analysis.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --time=48:00:00

# Load required modules
module purge
module load R


# Run the R script
R -f contact_analysis.R > contact_analysis.log 2> contact_analysis.err 