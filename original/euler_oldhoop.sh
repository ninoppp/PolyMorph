#!/bin/bash
#SBATCH --job-name=Oldhoop      # Job name    (default: sbatch)
#SBATCH --output=Oldhoop-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Oldhoop-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=128        # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=00:45:00           # Wall clock time limit

module load gcc
module list

echo "compiling old version ... "
g++ -fopenmp -O3 -o polyhoop_old.out polyhoop_old.cpp
OMP_NUM_THREADS=128
export OMP_NUM_THREADS
echo "running old polyhoop in parallel with $OMP_NUM_THREADS threads... "
./polyhoop_old.out

echo "all done."