#!/bin/bash
#SBATCH --job-name=genVarwidth      # Job name    (default: sbatch)
#SBATCH --output=genVarwidth-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=genVarwidth-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --ntasks=1
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=128       
#SBATCH --mem-per-cpu=1024        
#SBATCH --time=05:00:00         

module load gcc
module list

make clean
make generate_varwidth

export OMP_NUM_THREADS=128

echo "running patterning precision in parallel with $OMP_NUM_THREADS threads... "
srun ./polymorph $SLURM_NODEID # NodeID is used to determine individual seeds

echo "all done."