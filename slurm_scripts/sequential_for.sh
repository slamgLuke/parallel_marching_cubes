#!/bin/bash
#SBATCH -J Sequential # job name
#SBATCH -p big-mem # partition name
#SBATCH --mem=64G
#SBATCH --cpus-per-task 48 # number of CPUs per task
#SBATCH --account=investigacion1
#SBATCH -o log_sequential_for.out # output file

./sequential_for 32
echo "\n"
./sequential_for 64
echo "\n"
./sequential_for 128
echo "\n"
./sequential_for 256
echo "\n"
./sequential_for 512
echo "\n"
./sequential_for 1024
