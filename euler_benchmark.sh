#!/bin/bash
#SBATCH --job-name=benchmark_typical      # Job name    (default: sbatch)
#SBATCH --output=benchmark_typical-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=benchmark_typical-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1 
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=08:00:00           

module load gcc
module list

make clean
make benchmark_typical

# Run the program for OMP_NUM_THREADS equal to 1, 2, 4, 8, 16, 32, 64, 128
for ((i=4; i<=5; i++))
do
  OMP_NUM_THREADS=$((2**i))
  echo "Running benchmark with OMP_NUM_THREADS=$OMP_NUM_THREADS"
  export OMP_NUM_THREADS
  ./polymorph
done
echo "all done."