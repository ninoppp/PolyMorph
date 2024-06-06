#!/bin/bash
#SBATCH --job-name=Polymorph      # Job name    (default: sbatch)
#SBATCH --output=Polymorph-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Polymorph-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=16        # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=00:05:00           # Wall clock time limit

module load gcc
module list

echo "compiling polymorph ... "
make clean
make

export OMP_NUM_THREADS=16

echo "running polymorph in parallel with $OMP_NUM_THREADS threads... "
./polymorph

echo "all done."