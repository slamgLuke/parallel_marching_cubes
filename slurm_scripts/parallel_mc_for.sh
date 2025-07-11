#!/bin/bash
#SBATCH -J Parallel_for # job name
#SBATCH -p big-mem # partition name
#SBATCH --mem=64G
#SBATCH --cpus-per-task 48 # number of CPUs per task
#SBATCH --account=investigacion1
#SBATCH -o log_parallel_mc_for.out # output file

./parallel_mc_for 32
./parallel_mc_for 64
./parallel_mc_for 128
./parallel_mc_for 256
./parallel_mc_for 512
./parallel_mc_for 1024