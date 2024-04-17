#!/bin/bash
#SBATCH --job-name=<descriptivename>
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --time=0-00:30:00
#SBATCH --qos=short
#SBATCH --account=ageorge

module purge
module load cuda/10.1

make cuda

./build/boids_cuda