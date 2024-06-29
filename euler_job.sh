#!/bin/bash
#SBATCH --job-name=vel     # Job name    (default: sbatch)
#SBATCH --output=vel-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=vel-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=32       
#SBATCH --mem-per-cpu=2048        
#SBATCH --time=00:30:00     

module load gcc
module list

make clean
make velocity

export OMP_NUM_THREADS=32

echo "running job in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph

echo "all done."