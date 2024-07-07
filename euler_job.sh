#!/bin/bash
#SBATCH --job-name=typ_rect     # Job name    (default: sbatch)
#SBATCH --output=typ_rect-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=typ_rect-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=10
#SBATCH --nodes=10               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=8       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=08:00:00     

module load stack
module load gcc
module list

make clean
make benchmark_typical_rect

export OMP_NUM_THREADS=8

echo "running job in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph

echo "all done."