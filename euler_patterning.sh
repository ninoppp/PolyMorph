#!/bin/bash
#SBATCH --job-name=Patterning      # Job name    (default: sbatch)
#SBATCH --output=Patterning-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Patterning-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=120       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=03:00:00         

module load gcc
module list

echo "compiling patterning ... "
rm -f a.out
g++ -fopenmp -O3 -o a.out src/patterning_precision.cpp
OMP_NUM_THREADS=100
export OMP_NUM_THREADS
echo "running patterning precision in parallel with $OMP_NUM_THREADS threads... "
srun --exclusive ./a.out $SLURM_NODEID

echo "all done."