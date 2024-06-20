#!/bin/bash
#SBATCH --job-name=Benchmark      # Job name    (default: sbatch)
#SBATCH --output=Benchmark-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Benchmark-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=5
#SBATCH --nodes=5               
#SBATCH --ntasks-per-node=1 
#SBATCH --constraint=EPYC_7763    # Select node with CPU
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=12:00:00           

module load gcc
module list

make clean
make benchmark

# Run the program for OMP_NUM_THREADS equal to 1, 2, 4, 8, ..., 64, 128
for ((i=0; i<=7; i++))
do
  OMP_NUM_THREADS=$((2**i))
  echo "Running benchmark with OMP_NUM_THREADS=$OMP_NUM_THREADS"
  export OMP_NUM_THREADS
  ./polyhoop
done
echo "all done."