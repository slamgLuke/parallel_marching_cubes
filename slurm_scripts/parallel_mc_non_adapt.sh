#!/bin/bash
#SBATCH -J Parallel_non_adapt # job name
#SBATCH -p big-mem # partition name
#SBATCH --mem=64G
#SBATCH --cpus-per-task 48 # number of CPUs per task
#SBATCH --account=investigacion1
#SBATCH -o log_parallel_mc_non_adapt.out # output file

./parallel_mc_non_adapt 32
echo "\n"
./parallel_mc_non_adapt 64
echo "\n"
./parallel_mc_non_adapt 128
echo "\n"
./parallel_mc_non_adapt 256
echo "\n"
./parallel_mc_non_adapt 512
echo "\n"
./parallel_mc_non_adapt 1024