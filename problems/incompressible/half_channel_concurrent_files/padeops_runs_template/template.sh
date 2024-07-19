#!/bin/bash
#SBATCH --job-name=wide
#SBATCH -N 2
#SBATCH -n 64
#SBATCH --time=00:59:00
#SBATCH --partition=cpu
#SBATCH --exclude=yellowstone-cpu-3-1,yellowstone-cpu-11-20
#SBATCH --mail-type=BEGIN,END            # Mail when job starts and ends
#SBATCH --mail-user=youngin@mit.edu      # Email recipient


export inputFile="input_concurrent.dat"
export problemDir="/home/ctrsp-2024/youngin/mypadeops/PadeOps/problems/incompressible/half_channel_concurrent_files/padeops_runs"

cd /home/ctrsp-2024/youngin/mypadeops/PadeOps/build_opti/problems/incompressible
date
pwd

mpirun ./half_channel_concurrent $problemDir/$inputFile
