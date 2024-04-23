#!/bin/bash
#SBATCH --job-name=Benchmark      # Job name    (default: sbatch)
#SBATCH --output=Benchmark-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Benchmark-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=128         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=02:00:00           # Wall clock time limit

module load gcc
module list

echo "Running benchmark job... "
g++ -fopenmp -O3 -o benchmark.out benchmark.cpp

# Run the program for OMP_NUM_THREADS equal to 1, 2, 4, 8, ..., 64, 128
for ((i=0; i<=7; i++))
do
  OMP_NUM_THREADS=$((2**i))
  echo "Running with OMP_NUM_THREADS=$OMP_NUM_THREADS"
  export OMP_NUM_THREADS
  ./benchmark.out
done
echo "all done."