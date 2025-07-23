#!/bin/bash
#SBATCH -J bask                      #Job name
#SBATCH -p normal_q                  #Partition requested
#SBATCH --nodes=1                    #Request exclusive access to all resources on 1 node
#SBATCH --ntasks-per-node=128
#SBATCH --threads-per-core=1         #1 threads per core每核一个线程
#SBATCH --time=WALLTIME              #Time limit
#SBATCH --account=are_diagnost
#SBATCH --output=stdout              #file for batch script's standard output
#SBATCH --exclusive


LEVEL=13


module reset                         #Remove all modules from your environment
module load foss/2023a                      #Add a module to your environment

module load ImageMagick/7.1.1-15-GCCcore-12.3.0
module load FFmpeg/6.0-GCCcore-12.3.0

NAME=drop

srun --mpi=pmix -K1 -n $SLURM_NTASKS ./$NAME -m WALLTIME $LEVEL\  
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS

