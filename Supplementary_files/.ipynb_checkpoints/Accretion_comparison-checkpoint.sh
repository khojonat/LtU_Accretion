#!/bin/bash
#SBATCH --job-name=FFnoAGN
#SBATCH --output=logs/FFnoAGN_%j.out      # Output log file (%j = job ID)
#SBATCH --error=logs/FFnoAGN_%j.err       # Error log file
#SBATCH --time=1:00:00                 # Max runtime (HH:MM:SS)
#SBATCH --mail-type=ALL                 # Send email at begin and end
#SBATCH --mail-user=yja6qa@virginia.edu # Your email
#SBATCH --partition=standard            #parallel                       #standard
#SBATCH -A torrey-group

###SBATCH --export=NONE
#SBATCH --ntasks=3 		            # Max is 96
#SBATCH --mem-per-cpu=7500
###SBATCH --nodes=1  		                #--> Need parallel partition, replace --ntasks with these two
###SBATCH --ntasks-per-node=64 

module purge
source ~/load_arepo.sh
module load miniforge/24.3.0-py3.11
conda activate kho_env  #py3forge

# Run script
python Accretion_comparison.py
