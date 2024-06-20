#!/bin/bash
#SBATCH --job-name=gen_area      # Job name    (default: sbatch)
#SBATCH --output=gen_area-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=gen_area-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=8       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=01:00:00    
#SBATCH --constraint=EPYC_7763    # specify CPU for benchmarking

module load gcc
module list

make clean
make benchmark

export OMP_NUM_THREADS=8

echo "running job in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph

echo "all done."