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
echo "running polyhoop in parallel with 128 threads... "
OMP_NUM_THREADS=128 ./polyhoop.out

echo "moving output files..."
local timestamp=$(date +%Y-%m-%d_%H-%M)
mkdir $$timestamp
mv *.vtp $timestamp/
mv *.vts $timestamp/ 

echo "all done."