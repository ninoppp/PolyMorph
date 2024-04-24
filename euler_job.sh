#!/bin/bash
#SBATCH --job-name=Polymorph      # Job name    (default: sbatch)
#SBATCH --output=Polymorph-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Polymorph-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=4         # Number of CPUs per task
#SBATCH --mem-per-cpu=2048        # Memory per CPU
#SBATCH --time=08:00:00           # Wall clock time limit

module load gcc
module list

echo "compiling polymorph ... "
rm -f polymorph.out
g++ -fopenmp -O3 -o polymorph.out main.cpp
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
echo "running polymorph in parallel with $OMP_NUM_THREADS threads... "
./polymorph.out

echo "all done."