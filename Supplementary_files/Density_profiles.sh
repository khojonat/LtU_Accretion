#!/bin/bash
#SBATCH -J rho_profiles
#SBATCH -A torrey-group
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH -p standard
#SBATCH --output=rho_profiles_%j.out
#SBATCH --error=rho_profiles_%j.err
#SBATCH --mail-user=yja6qa@virginia.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniforge
conda activate kho_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python Density_profiles.py

