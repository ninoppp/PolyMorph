#!/bin/bash
#SBATCH --job-name=Polymorph      # Job name    (default: sbatch)
#SBATCH --output=Polymorph-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Polymorph-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=128        # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=00:20:00           # Wall clock time limit

module load gcc
module list

echo "compiling polyhoop ... "
g++ -fopenmp -O3 -o polyhoop.out polyhoop.cpp
OMP_NUM_THREADS=128
export OMP_NUM_THREADS
echo "running polyhoop in parallel with $OMP_NUM_THREADS threads... "
./polyhoop.out

echo "all done."