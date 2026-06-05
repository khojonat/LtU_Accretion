#!/bin/bash
#SBATCH -J lenient # Plot_smooth
#SBATCH -A torrey-group
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH -p standard
#SBATCH --output=logs/lenient_%j.out # Smooth_dist_%j.out
#SBATCH --error=logs/lenient_%j.err # Smooth_dist_%j.err
#SBATCH --mail-user=yja6qa@virginia.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniforge
conda activate kho_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python extract_lenient_data.py # Plot_dist_smooth.py
