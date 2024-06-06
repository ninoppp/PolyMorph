#!/bin/bash
#SBATCH --job-name=Patterning      # Job name    (default: sbatch)
#SBATCH --output=Patterning-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=Patterning-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=10
#SBATCH --nodes=10               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=128       
#SBATCH --mem-per-cpu=2048        
#SBATCH --time=05:00:00         

module load gcc
module list

echo "compiling patterning ... "
make clean
make generate_tissues

export OMP_NUM_THREADS=128

echo "running patterning precision in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph $SLURM_NODEID # NodeID is used to determine individual seeds

echo "all done."