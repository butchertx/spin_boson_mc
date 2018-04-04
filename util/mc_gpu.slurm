#!/bin/bash
#SBATCH --job-name=MPI_test
#SBATCH --partition=commons
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=3G
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --export=ALL

srun ~/class_mc/spin_boson_mpi.exe full_input.txt .
