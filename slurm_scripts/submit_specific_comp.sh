#!/bin/bash
#SBATCH -J Specific # Printing # 
#SBATCH -A torrey-group
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH -p standard
#SBATCH --output=specific_%j.out # printing.out   # 
#SBATCH --error=specific_%j.err # printing.err    # 
#SBATCH --mail-user=yja6qa@virginia.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniforge
conda activate kho_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# python printing_QM_RM.py
python Specific_components.py 
