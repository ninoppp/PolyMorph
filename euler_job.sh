#!/bin/bash
#SBATCH --job-name=poserr_cv      # Job name    (default: sbatch)
#SBATCH --output=poserr_cv-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=poserr_cv-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=10
#SBATCH --nodes=10               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=120       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=04:00:00         

module load gcc
module list

make clean
make poserr_cv

export OMP_NUM_THREADS=120

echo "running job in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph

echo "all done."