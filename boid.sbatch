#!/bin/bash
#SBATCH --job-name=<descriptivename>
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --cluster=smp
#SBATCH --time=0-00:30:00
#SBATCH --qos=short
#SBATCH --account=ageorge

module purge
module load gcc/12.2.0

## To add OpenMP support, you pass to the GNU compiler the -fopenmp option.
## The -o <desired_executable_name> instructs the compiler to name the binary 
## instead using the default a.out

make boidgen
./build/boids