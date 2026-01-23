#!/bin/bash
#SBATCH -J load_temp
#SBATCH -A torrey-group
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=7500
#SBATCH --time=24:00:00
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
