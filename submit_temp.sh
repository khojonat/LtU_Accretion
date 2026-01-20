#!/bin/bash
#SBATCH -J load_temp
#SBATCH -A torrey-group
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=10:00:00
#SBATCH -p standard
#SBATCH --output=logs/temps_%j.out
#SBATCH --error=logs/temps_%j.err
#SBATCH --mail-user=yja6qa@virginia.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniforge
conda activate kho_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python Gas_Temps.py
