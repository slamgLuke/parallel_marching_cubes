#!/bin/bash
#SBATCH -J Get_OMP_threads # job name
#SBATCH -p big-mem # partition name
#SBATCH --mem=64G
#SBATCH --cpus-per-task 48 # number of CPUs per task
#SBATCH --account=investigacion1
#SBATCH -o log_get_omp_threads.out # output file