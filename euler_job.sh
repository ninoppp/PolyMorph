#!/bin/bash
#SBATCH --job-name=bench_rect      # Job name    (default: sbatch)
#SBATCH --output=bench_rect-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=bench_rect-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=128       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=03:00:00     
#SBATCH --constraint=EPYC_7763    # Select node with CPU

module load gcc
module list

make clean
make benchmark_rect

export OMP_NUM_THREADS=128

echo "running job in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph

echo "all done."